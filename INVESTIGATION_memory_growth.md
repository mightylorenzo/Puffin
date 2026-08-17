# Investigation: per-rank memory grows with step count

Puffin's peak resident memory per MPI rank grows with the number of integration
steps, on a fixed mesh, with a fixed number of output dumps. Find out why, and
say whether it is a bug or expected behaviour.

**Status: resolved 2026-08-12.** Two independent components, one a genuine bug
(fixed) and one expected behaviour. See [Findings](#findings) for the answer; the
sections after it are the original brief, kept as a record of what was measured
and ruled out beforehand.

The brief was written on 2026-08-11 against `dev-test`; the investigation ran the
next day. Working tree
`/Users/lawc/gitroot/github.mightylorenzo.puffin/puffin/src`. Puffin is a Fortran
MPI FEL simulation code; build with `make -j8 -C build/`.

## Findings

### Component 1 — leaked MPI request handles (a genuine bug, fixed)

Three places in `puffin/lib/parallel/para_field.f90` call `mpi_issend` into a
single scalar `req` that is overwritten before anything waits on it. A request
that is never completed leaves the request object alive inside Open MPI for the
rest of the run, so every call leaks one or more:

| routine | bug | leaks on |
|---|---|---|
| `upd8a` | two `issend` blocks (`ac_rl`, `ac_il`) share one `mpi_wait` | every rank != 0, per call |
| `upd8da` periodic block | the size-handshake `issend` is overwritten by the next `issend` | rank `size-1`, per call |
| `calcBuff` | a loop over peer ranks issues `size-2` sends into one scalar handle | rank 0, per redistribution |

Fixed in commit `2035f4d` on `dev-test`, using one request per outstanding
`issend` plus `mpi_waitall`. The same commit stops `calcBuff` passing the literal
`0` to a non-blocking send — a literal may live in a compiler temporary that does
not outlive the call.

How it was found, in case the next one needs the same treatment:

1. **Sample RSS per rank during a single run** (15000 steps, 2 ranks, `ps` every
   0.3 s). Rank 1 climbed 35 -> 185 MB, linearly, with no plateau; rank 0 stayed
   flat at ~45 MB. Asymmetric growth kills the active-section hypothesis on its
   own — that would be bounded and would hit both ranks.
2. **`leaks <pid>`** reported only 25 KB unreachable but 80 MB of *live* malloc,
   so whatever was growing was still referenced. Not a classic leak.
3. **Two `heap <pid>` snapshots, diffed by class.** The entire delta was
   `malloc in opal_free_list_grow_st` — Open MPI's internal free list — at
   +32 MB / +513 blocks per 5 s, while every Puffin allocation stayed flat to
   within 0.3 MB. That is the accumulating request objects.

Both tools need `MallocStackLogging=1` in the environment to attribute blocks.

### Component 2 — the active section tracking the beam (expected, not a leak)

With the leak fixed, a heap diff on the big config shows the remaining growth
entirely in Puffin's own arrays:

```
+147.66 MB   +22 blocks   getLocalFieldIndices    (ac_/fr_/bk_ r+i field arrays)
+106.02 MB    +3 blocks   redist2FFTWlt           (tre_fft, tim_fft)
 +20.27 MB   +41 blocks   allact_rk4_arrs         (RK4 stage arrays)
```

Byte counts balloon while block counts stay tiny, so these arrays are being
*resized*, not accumulated. That is the third hypothesis below, confirmed: the
active section grows as the beam spreads in z2. `opal_free_list_grow_st` is now
static across the same interval (1.677 MB, 20 -> 21 blocks), confirming the leak
is gone on the big mesh too.

### Numbers

Small 3D config (`clara_test.in`), 2 ranks, peak RSS per rank. Component 1
dominates here, so the fix is decisive:

| steps | before | after |
|-------|--------|-------|
| 1500 | 43 MB | 31 MB |
| 15000 | 184 MB | 37 MB |
| 30000 | — | 39 MB |

Growth falls from ~10 kB/step to ~0.07 kB/step.

Big config (`shcse_big.in`), 6 ranks, dump cadences pinned to the step count.
Component 2 dominates late in the run, so the fix helps much less:

| steps | before | after |
|-------|--------|-------|
| 300 | 743 MB | 746 MB |
| 1200 | 785 MB | 751 MB |
| 3000 | 921 MB | 875 MB |

The leak was the constant-rate part: the early slope drops 47 -> 5.6 kB/step,
while the late slope barely moves (76 -> 69 kB/step) because that is component 2.

At 7500 steps with the big test's own cadence (`iWriteIntNthSteps = 750`), peak
RSS was 1493 MB before the fix and 1518 MB after — the difference is run-to-run
noise, because by then component 2 swamps the leak.

### Consequences for the big regression test

**~1.5 GB/rank at 7500 steps, and still climbing** — not the ~1.4 GB that
`test/testMPIIntegration3DBig.pf` is sized against. That estimate needs revising
upward, independently of this fix.

### Still open

The ceiling for component 2 is not established. It is bounded in principle — the
active section cannot exceed the mesh — but nothing here says where that bound
sits, so the memory needed by a run longer than 7500 steps is still unknown. Two
cheap ways to pin it down:

- Print `tllen` per rank at each redistribution and do one 3000-step run (~4 min)
  to see how much headroom is left before the active section spans the mesh.
- The discriminating experiment described below: a cold or low-energy-spread beam
  at high step count, to confirm memory tracks beam z2 extent rather than steps.

## The observation

Identical input, identical 85 x 85 x 5780 field mesh, identical rank count, only
`nPeriods` varied (30 steps per period):

| steps | peak RSS per rank |
|-------|-------------------|
|   300 |  743 MB |
|  1200 |  785 MB |
|  3000 |  921 MB |

Same effect at a different decomposition: at 16 ranks, 300 steps gave 352 MB and
7500 steps gave 1000 MB per rank. So it is not specific to a rank count, and the
per-rank growth is roughly 47-90 kB/step, with the apparent slope increasing at
larger step counts (i.e. possibly superlinear — worth confirming, my three points
are too few to fit anything with confidence).

## Reproduce it cheaply — use this, not the big run

The effect shows just as clearly on a small, fast config. `test/inputs/3D/clara_test.in`
is 3D (so it exercises the same `outputH5Field3DSD` writer), reads a frozen beam
from `clara_electrons_0.h5` so it is deterministic, and runs in **under a second**
at its default 150 steps. Raising `nPeriods` to 500 gives 15000 steps in ~8
seconds and reproduces the growth at 4.3x. Iterate here; do not use the 9-minute
big run for this investigation.

A sweep over both axes on that config (2 ranks, peak RSS of the largest rank via
`/usr/bin/time -l`):

| experiment | varied | held fixed | peak RSS |
|---|---|---|---|
| A | integrated dumps 3 -> 502 | 1500 steps | 43 -> 44 MB (flat) |
| B | field dumps 3 -> 152 | 1500 steps | 43 -> 43 MB (flat) |
| C | steps 1500 -> 15000 | 3 dumps | 43 -> 184 MB (4.3x) |

Baseline RSS with no integration to speak of is ~33 MB, so the growth in C is
large compared to the working set. Growth is ~7-10 kB/step and the slope rises
with step count (6.7, 8.3, 10.1 kB/step over the three points), i.e. mildly
superlinear — matching the big run, where it was 47-90 kB/step on a mesh ~1000x
larger.

The sweep script was ~30 lines of sed over `nPeriods`, `iWriteNthSteps` and
`iWriteIntNthSteps`. It lived in a session scratchpad under `/private/tmp` and is
gone; rewriting it takes a couple of minutes.

## Reproduction on the big config (for confirmation only)

The original script also lived in a `/private/tmp` scratchpad and is gone; it is
a few lines. Take
`test/inputs/3D/shcse_big.in`, copy it plus `test/inputs/3D/shcse_big_beam.in`
and `test/inputs/3D/seed_file.in` into `<dir>/inputs/3D/`, set `nPeriods` to the
value you want and set BOTH `iWriteNthSteps` and `iWriteIntNthSteps` equal to the
total step count (so the dump count stays fixed as you vary steps), then:

```
cd <dir>
OMP_NUM_THREADS=1 /usr/bin/time -l mpiexec -n 6 <src>/build/puffin/puffin inputs/3D/shcse_big.in
```

`maximum resident set size` in the `time -l` output is the largest single rank.

Practical notes, learned the hard way:
- `OMP_NUM_THREADS=1` always. OpenMP+MPI oversubscribes badly on Apple Silicon
  (up to 7x slower and very noisy).
- MPI runs need the Bash tool's `dangerouslyDisableSandbox: true` or they
  spuriously segfault.
- Each run writes ~1.3 GB (two 637 MB field dumps). Clean up between runs.
- Wall clock at 6 ranks: ~90 s for 300 steps, ~250 s for 3000, ~8.5 min for 7500.
  This machine has 6 performance + 12 efficiency cores and 48 GB.

## What the measurement does and does not show

**There is a growth component tied to step count.** All three runs above wrote
exactly 9 files (3 `aperp`, 3 `electrons`, 3 `integrated`), because the write
cadences were pinned to the total step count. Dump count identical across all
three, memory still growing — so something scales with steps, independently of
how often output is written.

**A per-dump leak was looked for directly and not found.** Experiments A and B
above vary the dump count with the step count held fixed — the mirror image of
the step sweep, and the test that would expose per-dump leakage. Peak RSS is flat
from 3 to 502 integrated dumps and from 3 to 152 field dumps. Whatever is growing
is not driven by the number of output writes, in this configuration.

This matters because the prior suspicion was an HDF5 library issue that has
bitten this codebase before, observed as memory climbing with each dump. One
reconciliation: if memory climbs steadily with *steps*, and dumps are simply when
one looks at it, steady growth would present as growth "at each dump" — the dumps
being the sampling points rather than the cause. The earlier sighting may also
predate this code, or involve the 1D writer or the `qDump` restart path, neither
of which the sweep above exercised. If you want to close that off, repeat
experiment A against those paths specifically.

If the HDF5 angle is pursued anyway: `h5open_f`/`h5close_f` are called per file
read/write throughout `puffin/lib/io/`, and HDF5 object IDs that are not
explicitly closed accumulate. `H5Fget_obj_count` on the file ID before close is a
quick way to see leaked IDs.

## What is ruled out

**It is not macroparticle imbalance across ranks.** The natural guess is that the
beam occupies a narrow z2 band so one rank ends up holding it all. It does not
work that way: Puffin defines front, back and active sections of the field mesh,
where the active section is where the beam is, and each of those three sections
is separately parallelised into slabs. (Confirmed by the code owner.)

## A third hypothesis for the step-linked component

*(Confirmed — this is component 2 above.)*

The active section tracks the beam, and the beam spreads and slips in z2 over the
run. So the active section grows with time while front/back shrink. If the
arrays or buffers for the three sections are grown but never shrunk — or if a
buffer is reallocated to a high-water mark — total allocation would climb with
step count in exactly this way, and would be a design consequence rather than a
leak. Relevant code:

- `puffin/lib/parallel/para_field.f90` — `getLocalFieldIndices`, `calcBuff`
- `puffin/lib/interaction/` — the RK4 chain (`puffin_mpi_RK4.f90`, `derivative.f90`, `rhs.f90`)
- redistribution is controlled by `sRedistLen` / `iRedistStp` (not set in
  `shcse_big.in`, so defaults apply — worth checking what those defaults do.
  Note `clara_test.in` *does* set them, `sRedistLen = 0.25` and
  `iRedistStp = 10`, so the two configs do not redistribute on the same cadence)

Distinguishing the two: if it is the active-section growth, memory should track
the *beam's z2 extent*, not the step count directly — so a run with a cold,
non-spreading beam (or `q_noise = .false.`, or a much smaller energy spread)
should flatten the curve while keeping the step count high. That is a clean
discriminating experiment and cheaper than instrumenting allocations.

## What a good answer looks like

*(Answered — see [Findings](#findings). It turned out to be both 1 and 2 at once:
a genuine leak that dominates small meshes and the early part of big runs, plus
expected active-section growth that dominates late on the big mesh. Peak RSS at
7500 steps is ~1.5 GB/rank and does not plateau; the ceiling is still unknown.)*

Identify the specific allocation(s) that grow, with evidence (a heap profiler, or
instrumented `allocated()`/size reporting around the suspects, or `leaks`/
`heap` on macOS). Then say whether it is:

1. a genuine leak (something allocated per step and never freed) — fix it; or
2. bounded, expected growth from the active section tracking the beam — in which
   case say what the actual ceiling is for a long run, since that sets the memory
   requirement for the big regression test.

Either way the practical question to answer is: what is peak RSS per rank for a
7500-step run, and does it plateau or keep climbing? The new
`test/testMPIIntegration3DBig.pf` is sized on the assumption that ~1.4 GB/rank at
7500 steps is the ceiling; if it keeps climbing, longer runs will need a revised
estimate.

## Environment

- Apple M5 Pro, 18 logical cores (6 performance + 12 efficiency), 48 GiB
- GNU Fortran (Homebrew GCC 16.1.0), Open MPI 5.0.9
- HDF5 from `/opt/homebrew/opt/hdf5-mpi` (parallel build), h5dump 2.2.0
- Build: `make -j8 -C build/`; tests: `ctest --test-dir build/`
