#!/usr/bin/env python3
"""Launch Puffin runs in the background and report on them.

One `Run` is one invocation: its own directory holding the four input files,
the log, and every HDF5 dump the run produces. Nothing is shared between
runs and nothing is overwritten, so a result can always be traced back to
the exact deck that produced it — the directory *is* the record.

Progress comes from counting `*_integrated_*.h5` dumps against the step
count the deck implies, not from parsing stdout. Puffin writes those on a
fixed cadence, so the count is a true fraction of work done; stdout only
announces module boundaries. Where the step count cannot be determined the
run reports as indeterminate rather than inventing a denominator.

`OMP_NUM_THREADS=1` is forced unless the caller overrides it. Threading is a
net loss with the default libgomp build, and leaving it unset lets OpenMP
claim every core on top of the MPI ranks, which oversubscribes badly.
"""

import os
import re
import glob
import time
import signal
import shutil
import subprocess

STATUS_READY = 'ready'
STATUS_RUNNING = 'running'
STATUS_DONE = 'done'
STATUS_FAILED = 'failed'
STATUS_STOPPED = 'stopped'

# Puffin's own end-of-run verdicts, printed by app/main.f90.
_OK_LINE = 'Puffin simulation completed'
_FAIL_LINE = 'Puffin simulation failed'


def find_puffin(explicit=None):
    """Locate the puffin binary: argument, then $PUFFIN_BIN, then the build tree."""
    for cand in (explicit, os.environ.get('PUFFIN_BIN')):
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, '..', '..'))
    # The build layout varies with the generator; .exe covers a Cygwin or
    # MinGW build, where the rest of the POSIX story still holds.
    for rel in ('build/puffin/puffin', 'build/bin/puffin', '../build/puffin/puffin'):
        for cand in (os.path.join(repo, rel), os.path.join(repo, rel + '.exe')):
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return os.path.abspath(cand)
    found = shutil.which('puffin')
    return os.path.abspath(found) if found else None


def find_mpiexec():
    for name in ('mpiexec', 'mpirun'):
        found = shutil.which(name)
        if found:
            return found
    return None


def slugify(text, fallback='run'):
    slug = re.sub(r'[^A-Za-z0-9._-]+', '-', (text or '').strip()).strip('-')
    return slug[:48] or fallback


class Run:
    """A single Puffin invocation in its own directory."""

    def __init__(self, directory, main_file, ranks=1, puffin_bin=None,
                 expected_steps=None, write_int_every=None, env=None):
        self.dir = os.path.abspath(directory)
        self.main_file = main_file
        self.ranks = max(1, int(ranks))
        self.puffin_bin = puffin_bin
        self.expected_steps = expected_steps
        self.write_int_every = write_int_every
        self.env_extra = dict(env or {})

        self.proc = None
        self.status = STATUS_READY
        self.started_at = None
        self.finished_at = None
        self.returncode = None
        self.error = None
        self._log_path = os.path.join(self.dir, 'run.log')
        self._log_fh = None

    # ── derived facts ───────────────────────────────────────────────────

    @property
    def name(self):
        return os.path.basename(self.dir)

    @property
    def log_path(self):
        return self._log_path

    @property
    def expected_dumps(self):
        """Scheduled integrated dumps: one at step 0, then every cadence steps.

        Puffin may also write a final end-of-run dump one step past the last
        scheduled one (`qDumpEnd`), so a finished run can leave one *more*
        file than this. That is why `progress` clamps instead of modelling
        it: erring this way means the bar reaches full at or just before the
        end, rather than sticking at 97% on a run that has actually finished.
        """
        if not self.expected_steps or not self.write_int_every:
            return None
        return int(self.expected_steps // self.write_int_every) + 1

    def command(self):
        """The argv this run will execute."""
        binary = self.puffin_bin or find_puffin()
        if not binary:
            raise RuntimeError(
                'puffin binary not found — set PUFFIN_BIN or pass --args puffin=/path/to/puffin')
        if self.ranks > 1:
            mpi = find_mpiexec()
            if not mpi:
                raise RuntimeError('mpiexec not found, cannot run on %d ranks' % self.ranks)
            return [mpi, '-n', str(self.ranks), binary, self.main_file]
        return [binary, self.main_file]

    def command_text(self):
        try:
            argv = self.command()
        except RuntimeError as exc:
            return '# %s' % exc
        env = self.environment_overrides()
        prefix = ' '.join('%s=%s' % kv for kv in sorted(env.items()))
        return ('%s %s' % (prefix, ' '.join(argv))).strip()

    def environment_overrides(self):
        env = {'OMP_NUM_THREADS': '1'}
        env.update(self.env_extra)
        return env

    # ── lifecycle ───────────────────────────────────────────────────────

    def start(self):
        if self.status == STATUS_RUNNING:
            raise RuntimeError('run already in progress')
        argv = self.command()
        env = dict(os.environ)
        env.update(self.environment_overrides())
        self._log_fh = open(self._log_path, 'w', buffering=1)
        self._log_fh.write('# %s\n# in %s\n\n' % (self.command_text(), self.dir))
        try:
            self.proc = subprocess.Popen(
                argv, cwd=self.dir, stdout=self._log_fh,
                stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as exc:
            self.status = STATUS_FAILED
            self.error = str(exc)
            self._close_log()
            raise
        self.status = STATUS_RUNNING
        self.started_at = time.time()
        return self

    def poll(self):
        """Refresh status from the process. Safe to call on any state."""
        if self.status != STATUS_RUNNING or self.proc is None:
            return self.status
        code = self.proc.poll()
        if code is None:
            return self.status
        self.returncode = code
        self.finished_at = time.time()
        self._close_log()
        if self.status == STATUS_STOPPED:
            return self.status
        # Puffin can exit 0 having logged a failure, so trust the log line
        # over the return code when the two disagree.
        verdict = self._log_verdict()
        if code == 0 and verdict is not False:
            self.status = STATUS_DONE
        else:
            self.status = STATUS_FAILED
            if verdict is False:
                self.error = 'Puffin reported failure — see run.log'
            elif code != 0:
                self.error = 'exited with code %d' % code
        return self.status

    def _log_verdict(self):
        """True if the log says completed, False if failed, None if neither."""
        tail = self.log_tail(4000)
        if _FAIL_LINE in tail:
            return False
        if _OK_LINE in tail:
            return True
        return None

    def stop(self):
        """Terminate the run and everything mpiexec spawned.

        Signalling the whole process group is what matters: killing mpiexec
        alone can leave the ranks running. `start_new_session=True` at launch
        put them in their own group so this cannot reach anything else.

        Windows has no process groups in this sense (and no killpg, getpgid
        or SIGKILL), so there we fall back to terminating mpiexec directly.
        Puffin on Windows is built under WSL or Cygwin, which are POSIX and
        take the group path.
        """
        if self.proc is None or self.proc.poll() is not None:
            return
        self.status = STATUS_STOPPED
        self._signal_group(signal.SIGTERM, self.proc.terminate)
        try:
            self.proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self._signal_group(getattr(signal, 'SIGKILL', signal.SIGTERM),
                               self.proc.kill)
        self.finished_at = time.time()
        self._close_log()

    def _signal_group(self, sig, fallback):
        try:
            os.killpg(os.getpgid(self.proc.pid), sig)
        except AttributeError:          # no process groups (Windows)
            fallback()
        except (ProcessLookupError, PermissionError):
            fallback()

    def _close_log(self):
        if self._log_fh and not self._log_fh.closed:
            self._log_fh.close()
        self._log_fh = None

    # ── observation ─────────────────────────────────────────────────────

    def dump_count(self):
        return len(glob.glob(os.path.join(self.dir, '*_integrated_*.h5')))

    def progress(self):
        """Completed fraction in 0..1, or None when it cannot be known."""
        total = self.expected_dumps
        if not total:
            return None
        return min(1.0, self.dump_count() / float(total))

    def elapsed(self):
        if self.started_at is None:
            return 0.0
        end = self.finished_at if self.finished_at else time.time()
        return end - self.started_at

    def log_tail(self, nbytes=8000):
        try:
            size = os.path.getsize(self._log_path)
            with open(self._log_path, 'r', errors='replace') as fh:
                if size > nbytes:
                    fh.seek(size - nbytes)
                    fh.readline()          # drop the partial first line
                return fh.read()
        except OSError:
            return ''

    def current_module(self):
        """Which undulator module Puffin last announced, if any."""
        matches = re.findall(r'Simulating undulator module\s+(\d+)', self.log_tail())
        return int(matches[-1]) if matches else None

    def field_dumps(self):
        return sorted(glob.glob(os.path.join(self.dir, '*_aperp_*.h5')))

    def is_1d(self, main_deck=None):
        """Whether to open the 1D or 3D viewer for this run's output."""
        if main_deck is not None:
            return bool(main_deck.get('qOneD', False))
        for path in self.field_dumps()[:1]:
            try:
                import h5py
                with h5py.File(path, 'r') as fh:
                    shape = fh['aperp'].shape
                    return len(shape) <= 2 or (shape[-1] == 1 and shape[-2] == 1)
            except Exception:
                return False
        return False
