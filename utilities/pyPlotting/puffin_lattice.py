#!/usr/bin/env python3
"""Model of a Puffin lattice file, for reading, editing and writing back.

A lattice is an ordered sequence of elements, each a tag and a row of
positional numbers:

    QU  fx fy                                        quadrupole
    UN  type periods alpha taper steps ux uy kbnx kbny   undulator
    CH  zbar slip disp                               chicane
    DR  zbar                                         drift
    MO  wavenum mag                                  energy modulation

The order *is* the machine. An element has no position of its own: its z is
the sum of the lengths before it, so the sequence is the only thing an
editor can meaningfully rearrange. Reordering and resizing are edits;
"moving an element to z = 3 m" is not a thing the format can express.

Comments are kept. Those above the first element are the file header and
stay at the top; comments directly above an element belong to it and travel
with it when it moves, which is what makes a reordered file still read like
the one you started with.
"""

import re

# tag → (field names, human label). Field order is the read order in
# acc_lattice.f90; getting it wrong silently changes the machine.
SPEC = {
    'UN': (['type', 'periods', 'alpha', 'taper', 'steps', 'ux', 'uy', 'kbnx', 'kbny'],
           'Undulator'),
    'DR': (['zbar'], 'Drift'),
    'QU': (['fx', 'fy'], 'Quad'),
    'CH': (['zbar', 'slip', 'disp'], 'Chicane'),
    'MO': (['wavenum', 'mag'], 'Modulation'),
}

ORDER = ['UN', 'DR', 'QU', 'CH', 'MO']

# Sensible starting values for a newly added element.
DEFAULTS = {
    'UN': ["'planepole'", '10', '1.0', '0.0', '30', '1.0', '1.0', '0.0', '0.0'],
    'DR': ['1.0'],
    'QU': ['1.0', '-1.0'],
    'CH': ['0.0', '0.0', '0.0'],
    'MO': ['1.0', '0.0'],
}

_NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')


def _tokenize(text):
    """Split a lattice line into tokens, keeping quoted strings whole."""
    return re.findall(r"'[^']*'|\"[^\"]*\"|\S+", text)


class Element:
    """One lattice line: a tag, its positional fields, and its comments."""

    def __init__(self, tag, fields, comment='', above=None):
        self.tag = tag
        self.fields = list(fields)
        self.comment = comment          # trailing `! ...` on the line
        self.above = list(above or [])  # comment lines owned by this element

    # ── field access ────────────────────────────────────────────────────

    @property
    def names(self):
        return SPEC.get(self.tag, ([], ''))[0]

    @property
    def label(self):
        return SPEC.get(self.tag, ([], self.tag))[1]

    def get(self, name, default=None):
        try:
            return self.fields[self.names.index(name)]
        except (ValueError, IndexError):
            return default

    def set(self, name, value):
        try:
            i = self.names.index(name)
        except ValueError:
            return False
        while len(self.fields) <= i:
            self.fields.append('0.0')
        self.fields[i] = str(value).strip()
        return True

    def number(self, name, default=0.0):
        """A field as a float, or `default` if it is absent or not numeric."""
        raw = self.get(name)
        if raw is None:
            return default
        try:
            return float(str(raw).replace('d', 'e').replace('D', 'e'))
        except ValueError:
            return default

    # ── geometry ────────────────────────────────────────────────────────

    def length(self, lambda_w):
        """Displayed length in metres, matching the viewers' convention.

        The viewers scale periods and zbar by lambda_w to lay the diagram
        out; this follows them exactly so the editor and the result panels
        draw the same machine the same way.
        """
        if self.tag == 'UN':
            return self.number('periods') * lambda_w
        if self.tag in ('DR', 'CH'):
            return self.number('zbar') * lambda_w
        return 0.0

    def steps(self):
        """Integration steps this element contributes (undulators only)."""
        if self.tag != 'UN':
            return 0
        return int(self.number('periods')) * int(self.number('steps'))

    # ── text ────────────────────────────────────────────────────────────

    def line(self):
        body = '%-3s %s' % (self.tag, '  '.join(str(f) for f in self.fields))
        return (body + '  ' + self.comment).rstrip() if self.comment else body

    def copy(self):
        return Element(self.tag, self.fields, self.comment, self.above)

    def __repr__(self):
        return '<%s %s>' % (self.tag, ' '.join(str(f) for f in self.fields))


class Lattice:
    """An ordered list of elements, plus the file header."""

    def __init__(self, elements=None, header=None):
        self.elements = list(elements or [])
        self.header = list(header or [])

    # ── parsing ─────────────────────────────────────────────────────────

    @classmethod
    def parse(cls, text):
        header, pending, elements = [], [], []
        for raw in (text or '').splitlines():
            stripped = raw.strip()
            code = stripped.split('!')[0].strip()
            # Case-sensitive, exactly as acc_lattice.f90 compares it. Being
            # more permissive here would be worse than useless: an unmarked
            # header line like "Module file for Puffin" starts with "Mo",
            # which Fortran does not read as an element — matching it
            # case-insensitively would show a machine that is not the one
            # that runs. Two shipped decks depend on this.
            tag = code[:2] if code else ''
            if tag in SPEC:
                tokens = _tokenize(code)
                comment = stripped[len(stripped.split('!')[0]):].rstrip() \
                    if '!' in stripped else ''
                elements.append(Element(tag, tokens[1:], comment, pending))
                pending = []
            elif not elements and not pending:
                # Still in the file header until something attaches below.
                (header if stripped or header else header).append(raw)
            else:
                pending.append(raw)
        lat = cls(elements, header)
        lat._trailing = pending
        return lat

    def text(self):
        out = list(self.header)
        if out and out[-1].strip():
            out.append('')
        for el in self.elements:
            out.extend(el.above)
            out.append(el.line())
        out.extend(getattr(self, '_trailing', []))
        return '\n'.join(out).rstrip('\n') + '\n'

    # ── geometry ────────────────────────────────────────────────────────

    def layout(self, lambda_w):
        """[(element, z_start, length)] with z accumulated along the sequence."""
        out, z = [], 0.0
        for el in self.elements:
            L = el.length(lambda_w)
            out.append((el, z, L))
            z += L
        return out

    def total_length(self, lambda_w):
        return sum(el.length(lambda_w) for el in self.elements)

    def total_steps(self):
        return sum(el.steps() for el in self.elements)

    # ── editing ─────────────────────────────────────────────────────────

    def move(self, index, new_index):
        """Move the element at `index` to `new_index`. Returns where it landed."""
        if not (0 <= index < len(self.elements)):
            return index
        new_index = max(0, min(len(self.elements) - 1, new_index))
        if new_index == index:
            return index
        el = self.elements.pop(index)
        self.elements.insert(new_index, el)
        return new_index

    def insert(self, tag, index=None):
        if tag not in SPEC:
            raise ValueError('unknown element %r' % tag)
        el = Element(tag, DEFAULTS[tag])
        index = len(self.elements) if index is None else max(0, min(len(self.elements), index))
        self.elements.insert(index, el)
        return index

    def duplicate(self, index):
        if not (0 <= index < len(self.elements)):
            return index
        self.elements.insert(index + 1, self.elements[index].copy())
        return index + 1

    def delete(self, index):
        if 0 <= index < len(self.elements):
            self.elements.pop(index)
        return max(0, min(index, len(self.elements) - 1))

    def index_at(self, z, lambda_w):
        """Which element covers position `z`, or the nearest one."""
        best, best_d = None, None
        for i, (el, zs, L) in enumerate(self.layout(lambda_w)):
            if L > 0 and zs <= z < zs + L:
                return i
            centre = zs + L / 2.0
            d = abs(z - centre)
            if best_d is None or d < best_d:
                best, best_d = i, d
        return best

    def insertion_index(self, z, lambda_w):
        """Where a drop at `z` should place an element in the sequence."""
        for i, (el, zs, L) in enumerate(self.layout(lambda_w)):
            if z < zs + (L / 2.0 if L > 0 else 0.0):
                return i
        return len(self.elements)

    def __len__(self):
        return len(self.elements)
