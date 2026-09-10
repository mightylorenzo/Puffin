#!/usr/bin/env python3
"""Read, edit and write Puffin input decks without disturbing them.

A deck is a Fortran namelist file: `&MDATA ... /` blocks of `key = value`
lines, with `!` comments carrying the parameter documentation. Those comment
headers are the only reference material for most parameters, so an editor
that dropped them would make the file worse than it found it.

So this is deliberately not a parse-to-dict-and-regenerate module. A file is
held as its original lines, plus an index of where each assignment lives.
Editing a value rewrites the value text on that one line and leaves
everything else — comments, blank lines, spacing, key order, keys this
module has never heard of — exactly as it was. Round-tripping a file you
have not edited returns it byte for byte.

Fortran namelist rules that matter here:

  * keys are case-insensitive, so they are indexed casefolded
  * a repeated key takes its last value (`clara.in` really does set
    `meshType` twice), so `set` rewrites the last occurrence
  * subscripted keys (`sSigmaE(1,:)`, `dist_f(2)`) are distinct entries and
    are kept verbatim rather than being folded into an array

Values are typed on the way out (`.true.` → True, `1e-6` → float, comma
lists → list) and re-formatted on the way in, preserving Fortran spelling
for logicals and quoted strings.
"""

import re

_ASSIGN = re.compile(r'^(\s*)([A-Za-z_][A-Za-z_0-9]*(?:\([^)]*\))?)(\s*=\s*)(.*)$')
_BLOCK_OPEN = re.compile(r'^\s*&\s*([A-Za-z_][A-Za-z_0-9]*)')
_BLOCK_CLOSE = re.compile(r'^\s*/\s*$')


def split_comment(text):
    """Split a namelist line's value from its trailing `!` comment.

    A `!` inside quotes is data, not a comment — `zundType = 'plane!pole'`
    would otherwise lose its tail.
    """
    quote = None
    for i, ch in enumerate(text):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == '!':
            return text[:i], text[i:]
    return text, ''


def parse_value(text):
    """Fortran namelist value text → Python value.

    Returns a list for comma-separated values, a scalar otherwise. Anything
    unrecognised comes back as the stripped string, so an odd value survives
    a load/save cycle rather than raising.
    """
    text = text.strip().rstrip(',').strip()
    if not text:
        return ''
    parts = _split_commas(text)
    vals = [_parse_scalar(p) for p in parts]
    return vals[0] if len(vals) == 1 else vals


def _split_commas(text):
    """Split on commas that are not inside quotes."""
    out, buf, quote = [], '', None
    for ch in text:
        if quote:
            buf += ch
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            buf += ch
        elif ch == ',':
            out.append(buf)
            buf = ''
        else:
            buf += ch
    if buf.strip():
        out.append(buf)
    return [p.strip() for p in out if p.strip()]


def _parse_scalar(text):
    text = text.strip()
    low = text.lower()
    if low in ('.true.', '.t.', 't', 'true'):
        return True
    if low in ('.false.', '.f.', 'f', 'false'):
        return False
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    # Fortran writes d-exponents (1.0d-6); Python needs an e.
    numeric = text.replace('D', 'e').replace('d', 'e')
    try:
        return int(numeric)
    except ValueError:
        pass
    try:
        return float(numeric)
    except ValueError:
        return text


def format_value(value):
    """Python value → Fortran namelist value text."""
    if isinstance(value, (list, tuple)):
        return ', '.join(format_value(v) for v in value)
    if isinstance(value, bool):
        return '.true.' if value else '.false.'
    if isinstance(value, str):
        # A bare string that already looks like a namelist literal (a number,
        # a logical, an already-quoted string) is passed through as typed;
        # anything else is a character value and gets quoted.
        stripped = value.strip()
        if not stripped:
            return "''"
        if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in "'\"":
            return stripped
        if stripped.lower() in ('.true.', '.false.'):
            return stripped
        try:
            float(stripped.replace('D', 'e').replace('d', 'e'))
            return stripped
        except ValueError:
            pass
        if ',' in stripped:          # already a list the caller formatted
            return stripped
        return "'%s'" % stripped
    return repr(value)


class Deck:
    """One namelist file, editable in place.

    `deck['srho']` reads, `deck['srho'] = 0.005` writes. Keys are
    case-insensitive. Setting a key that is not present appends it inside
    the preferred block (the first block, unless `block=` says otherwise).
    """

    def __init__(self, lines, path=None):
        self.lines = list(lines)
        self.path = path
        self._reindex()

    # ── construction ────────────────────────────────────────────────────

    @classmethod
    def load(cls, path):
        with open(path, 'r', errors='replace') as fh:
            return cls(fh.read().splitlines(), path=path)

    @classmethod
    def loads(cls, text, path=None):
        return cls(text.splitlines(), path=path)

    def _reindex(self):
        """Map casefolded key → [line indices], and record block extents."""
        self._index = {}
        self.blocks = {}          # casefolded name → (open idx, close idx)
        block, opened = None, None
        for i, line in enumerate(self.lines):
            code, _ = split_comment(line)
            m_open = _BLOCK_OPEN.match(code)
            if m_open:
                block, opened = m_open.group(1).casefold(), i
                continue
            if _BLOCK_CLOSE.match(code):
                if block is not None:
                    self.blocks[block] = (opened, i)
                block, opened = None, None
                continue
            m = _ASSIGN.match(code)
            if m:
                self._index.setdefault(m.group(2).casefold(), []).append(i)

    # ── reading ─────────────────────────────────────────────────────────

    def __contains__(self, key):
        return key.casefold() in self._index

    def keys(self):
        """Keys in file order, as originally spelled."""
        seen = []
        for i, line in enumerate(self.lines):
            code, _ = split_comment(line)
            m = _ASSIGN.match(code)
            if m:
                seen.append(m.group(2))
        return seen

    def __getitem__(self, key):
        idx = self._index.get(key.casefold())
        if not idx:
            raise KeyError(key)
        code, _ = split_comment(self.lines[idx[-1]])   # last wins
        return parse_value(_ASSIGN.match(code).group(4))

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def raw(self, key, default=None):
        """The value text as written, without parsing."""
        idx = self._index.get(key.casefold())
        if not idx:
            return default
        code, _ = split_comment(self.lines[idx[-1]])
        return _ASSIGN.match(code).group(4).strip().rstrip(',').strip()

    # ── writing ─────────────────────────────────────────────────────────

    def __setitem__(self, key, value):
        self.set(key, value)

    def set(self, key, value, block=None):
        """Set `key`, rewriting the last occurrence or appending to a block."""
        text = format_value(value)
        idx = self._index.get(key.casefold())
        if idx:
            i = idx[-1]
            code, comment = split_comment(self.lines[i])
            m = _ASSIGN.match(code)
            trail = ',' if m.group(4).rstrip().endswith(',') else ''
            self.lines[i] = '%s%s%s%s%s%s' % (
                m.group(1), m.group(2), m.group(3), text, trail, comment)
            return
        self._append(key, text, block)

    def _append(self, key, text, block=None):
        """Insert `key = text` just before a block's closing `/`."""
        if not self.blocks:
            raise ValueError('no namelist block to add %r to' % key)
        if block is not None:
            target = self.blocks.get(block.casefold())
            if target is None:
                raise KeyError('no &%s block in %s' % (block, self.path))
        else:
            target = min(self.blocks.values(), key=lambda b: b[0])
        self.lines.insert(target[1], ' %s = %s' % (key, text))
        self._reindex()

    def update(self, mapping, block=None):
        for key, value in mapping.items():
            self.set(key, value, block=block)

    # ── output ──────────────────────────────────────────────────────────

    def text(self):
        return '\n'.join(self.lines) + '\n'

    def save(self, path):
        with open(path, 'w') as fh:
            fh.write(self.text())
        return path

    def __repr__(self):
        return '<Deck %s: %d keys>' % (self.path or '<memory>', len(self._index))


# ── lattice files ───────────────────────────────────────────────────────────

# A lattice line is an element code plus positional fields, e.g.
#   UN 'planepole'  29  1.0  0.0  30  1.0  1.0  0.0  0.0
#      type      periods alpha taper steps/period ux uy kbnx kbny
# Only the undulators consume integration steps; drifts, quads, chicanes and
# modulations are applied as instantaneous maps between them.
_UN_PERIODS, _UN_STEPS = 2, 5


def lattice_steps(path):
    """Total integration steps in a lattice file, or None if unreadable.

    Sum over undulator elements of periods × steps-per-period. This is what
    turns a dump count into a progress fraction.
    """
    try:
        with open(path, 'r', errors='replace') as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    total = 0
    for line in lines:
        code = line.split('!')[0].strip()
        if not code[:2].upper() == 'UN':
            continue
        fields = code.replace(',', ' ').split()
        try:
            total += int(float(fields[_UN_PERIODS])) * int(float(fields[_UN_STEPS]))
        except (IndexError, ValueError):
            return None            # malformed: better no estimate than a wrong one
    return total or None


def total_steps(main_deck, lattice_path=None):
    """Integration steps a run will take, or None if it cannot be determined.

    With a lattice file the lattice decides; without one, `nPeriods` and
    `stepsPerPeriod` in the main deck do.
    """
    if lattice_path:
        steps = lattice_steps(lattice_path)
        if steps:
            return steps
        return None
    try:
        return int(main_deck['nPeriods']) * int(main_deck['stepsPerPeriod'])
    except (KeyError, TypeError, ValueError):
        return None
