#!/usr/bin/env python3
"""Visual lattice editor: drag to reorder, click to edit.

Draws the machine in the same language as the result viewers — undulators
blue, chicanes orange, quads aqua, drift the recessive chrome grey — and
makes it editable in place.

The strip is inline SVG in a Div rather than a Bokeh figure. A figure is
drawn to canvas, which means no rounded caps, no CSS hover state, no
transition when an element moves, and no crisp text inside a bar. SVG gives
all of those, at the cost of wiring the pointer handling by hand:

    DocumentReady  installs one delegated listener set, once
    composedPath   finds the bar under the pointer across the shadow root
                   boundary that Bokeh puts every widget behind
    a hidden       carries {action, index, ...} back as JSON; a sequence
    TextInput      number makes repeat events distinct values so the
                   Python-side on_change always fires

Two views, because one geometry cannot serve both jobs:

    To scale   element widths are their real z extent, so it reads as the
               machine and matches the viewers' lattice panel
    Sequence   every element the same width, so each one is big enough to
               grab

That distinction is not decoration. CLARA's lattice is 222 elements over
21 m: to scale, each of its 170 drifts is well under a pixel, and nothing
smaller than a pixel can be clicked, let alone dragged.

What dragging means here is fixed by the file format. An element has no
position of its own — z is the running sum of the lengths before it — so a
drag reorders the sequence and never sets a coordinate. The drop marker
shows the gap the element will land in.
"""

import json

from bokeh.models import (Div, Button, Select, TextInput, RadioButtonGroup,
                          TextAreaInput, CustomJS)
from bokeh.layouts import column, row
from bokeh.events import DocumentReady

import puffin_viz_theme as pvt
import puffin_ui_style as pus
from puffin_lattice import Lattice, SPEC, ORDER

# Element colours by role, from the shared theme, and the bar heights the
# viewers use — undulators tallest, drift lowest, so the machine's shape is
# legible from the silhouette alone.
HEIGHT = {'UN': 1.00, 'CH': 0.76, 'DR': 0.46, 'QU': 0.88, 'MO': 0.88}
COLOUR_ROLE = {'UN': 'field', 'CH': 'beam', 'QU': 'select',
               'MO': 'ink2', 'DR': 'drift'}

EDIT_FIELDS = {
    'UN': ['type', 'periods', 'steps', 'alpha', 'taper'],
    'DR': ['zbar'],
    'QU': ['fx', 'fy'],
    'CH': ['zbar', 'slip', 'disp'],
    'MO': ['wavenum', 'mag'],
}

FIELD_HELP = {
    'periods': 'Undulator periods — sets both the length and the step count',
    'steps': 'Integration steps per period',
    'alpha': 'a_w scaling for this module, relative to the global saw',
    'taper': 'd(a_w)/dz within this module',
    'zbar': 'Length in scaled zbar',
    'slip': 'Slippage added by the chicane',
    'disp': 'Dispersion parameter',
    'fx': 'Focusing strength, x plane',
    'fy': 'Focusing strength, y plane',
    'wavenum': 'Modulation wavenumber',
    'mag': 'Modulation magnitude',
    'type': 'Undulator field model for this module',
}

# An empty type is not a missing value: it selects the 1D undulator with no
# off-axis variation of a_w, and the 1D decks rely on it. It has to be an
# option in its own right — without it the dropdown would show a nearby
# value like 'planepole' for those elements, and writing that back would
# quietly change the physics.
UND_TYPES = [("''", '1D — no off-axis variation'),
             ("'planepole'", 'planepole'),
             ("'curved'", 'curved'),
             ("'helical'", 'helical')]

STRIP_H = 92          # SVG height in px
BAR_TOP = 10
BAR_H = 54
MIN_W = 3.0           # narrowest a bar may draw, in px — below this it is
                      # not a pointer target and the view is lying about
                      # what can be grabbed


# The pointer wiring. Installed once per browser session; every strip on the
# page is found by class, so a reload or a second editor needs no new code.
_WIRE_JS = """
if (window.__pfLatticeWired) { return; }
window.__pfLatticeWired = true;

let seq = 0;
const send = (msg) => { msg.n = ++seq; bus.value = JSON.stringify(msg); };

const barAt = (ev) => {
  for (const n of ev.composedPath()) {
    if (n instanceof Element && n.classList && n.classList.contains('pf-bar')) return n;
  }
  return null;
};
const stripAt = (ev) => {
  for (const n of ev.composedPath()) {
    if (n instanceof Element && n.classList && n.classList.contains('pf-strip')) return n;
  }
  return null;
};
// Fraction across the strip, so Python can map to an insertion point using
// the same layout it drew from.
const fracAt = (ev, strip) => {
  const r = strip.getBoundingClientRect();
  return Math.max(0, Math.min(1, (ev.clientX - r.left) / r.width));
};

let drag = null;

document.addEventListener('pointerdown', (ev) => {
  const bar = barAt(ev);
  if (!bar) return;
  const strip = stripAt(ev);
  drag = {i: +bar.dataset.i, moved: false, strip: strip};
  bar.classList.add('pf-dragging');
  send({action: 'select', i: drag.i});
  ev.preventDefault();
}, true);

document.addEventListener('pointermove', (ev) => {
  if (!drag) {
    const bar = barAt(ev);
    document.querySelectorAll('.pf-bar.pf-hover').forEach(b =>
      b.classList.remove('pf-hover'));
    if (bar) bar.classList.add('pf-hover');
    return;
  }
  drag.moved = true;
  const f = fracAt(ev, drag.strip);
  // Move the drop marker in the browser so it tracks the pointer at frame
  // rate; the authoritative index still comes from Python on release.
  const mark = drag.strip.querySelector('.pf-drop');
  if (mark) {
    const w = drag.strip.viewBox.baseVal.width || drag.strip.clientWidth;
    mark.setAttribute('x1', f * w); mark.setAttribute('x2', f * w);
    mark.style.opacity = '1';
  }
}, true);

document.addEventListener('pointerup', (ev) => {
  if (!drag) return;
  const d = drag; drag = null;
  document.querySelectorAll('.pf-bar.pf-dragging').forEach(b =>
    b.classList.remove('pf-dragging'));
  if (d.moved) send({action: 'drop', i: d.i, f: fracAt(ev, d.strip)});
}, true);

document.addEventListener('pointercancel', () => {
  drag = null;
  document.querySelectorAll('.pf-bar.pf-dragging').forEach(b =>
    b.classList.remove('pf-dragging'));
}, true);
"""


class LatticeEditor:
    """Owns the lattice strip, its controls and the text that backs them."""

    def __init__(self, theme, on_change=None, width=940, lambda_w=0.0275,
                 doc=None):
        self.T = theme
        self.on_change = on_change or (lambda text: None)
        self.width = width
        self.lambda_w = lambda_w
        self.lat = Lattice.parse('')
        self.selected = None
        self.enabled = True
        self._suppress = False
        self._slots_cache = []

        self._build(doc)

    # ── construction ────────────────────────────────────────────────────

    def _build(self, doc):
        T = self.T
        self.icss = pus.input_css(T)

        self.bus = TextInput(value='', visible=False)
        self.bus.on_change('value', self._on_bus)

        self.strip = Div(text='', sizing_mode='stretch_width',
                         stylesheets=[pus.fill_css()])

        self.view_mode = RadioButtonGroup(
            labels=['Sequence', 'To scale'], active=0, width=168,
            stylesheets=[pus.radio_css(T)])
        self.view_mode.on_change('active', lambda a, o, n: self.redraw())

        self.summary = Div(sizing_mode='stretch_width',
                           stylesheets=[pus.fill_css()])

        self.add_select = Select(
            options=[(t, SPEC[t][1]) for t in ORDER], value='UN',
            width=132, stylesheets=[self.icss])
        self.add_btn = Button(label='＋ Add', width=76,
                              stylesheets=[pus.button_css(T, 'plain')])
        self.dup_btn = Button(label='Duplicate', width=92,
                              stylesheets=[pus.button_css(T, 'plain')])
        self.del_btn = Button(label='Delete', width=74,
                              stylesheets=[pus.button_css(T, 'stop')])
        self.left_btn = Button(label='◀', width=40,
                               stylesheets=[pus.button_css(T, 'plain')])
        self.right_btn = Button(label='▶', width=40,
                                stylesheets=[pus.button_css(T, 'plain')])

        self.add_btn.on_click(self._add)
        self.dup_btn.on_click(self._duplicate)
        self.del_btn.on_click(self._delete)
        self.left_btn.on_click(lambda: self._nudge(-1))
        self.right_btn.on_click(lambda: self._nudge(+1))

        self.sel_title = Div(sizing_mode='stretch_width',
                             stylesheets=[pus.fill_css()])
        self.field_row = row(sizing_mode='stretch_width', spacing=10)
        self.field_widgets = {}

        self.text = TextAreaInput(value='', rows=11,
                                  sizing_mode='stretch_width',
                                  stylesheets=[self.icss])
        self.text.on_change('value', self._on_text_edit)

        self.layout = column(
            row(self.view_mode, self.summary, spacing=14,
                sizing_mode='stretch_width'),
            self.strip,
            row(self.add_select, self.add_btn, self.dup_btn, self.del_btn,
                Div(text='', width=10), self.left_btn, self.right_btn,
                spacing=6),
            self.sel_title,
            self.field_row,
            Div(text=pus.hint(T, 'File text — edits here rebuild the diagram above.'),
                sizing_mode='stretch_width', stylesheets=[pus.fill_css()]),
            self.text,
            self.bus,
            sizing_mode='stretch_width', spacing=8)

        if doc is not None:
            doc.js_on_event(DocumentReady,
                            CustomJS(args=dict(bus=self.bus), code=_WIRE_JS))

    # ── external interface ──────────────────────────────────────────────

    def set_text(self, text, lambda_w=None):
        if lambda_w:
            self.lambda_w = lambda_w
        self.lat = Lattice.parse(text or '')
        self.selected = None
        self._suppress = True
        self.text.value = text or ''
        self._suppress = False
        self.redraw()

    def get_text(self):
        return self.lat.text()

    def set_enabled(self, enabled, message=''):
        self.enabled = enabled
        for w in (self.add_btn, self.dup_btn, self.del_btn, self.left_btn,
                  self.right_btn, self.add_select, self.text, self.view_mode):
            w.disabled = not enabled
        if not enabled:
            # Clear the model and the text too, not just the drawing. A deck
            # with no lattice that still showed the previous deck's file
            # would be describing a machine this run does not have.
            self.lat = Lattice.parse('')
            self.selected = None
            self._suppress = True
            self.text.value = ''
            self._suppress = False
            self.strip.text = self._empty_strip(
                message or 'This deck has no lattice file.')
            self.summary.text = ''
            self.sel_title.text = ''
            self.field_row.children = []

    # ── geometry ────────────────────────────────────────────────────────

    def _slots(self):
        """[(x, w)] in px across the strip's own coordinate space."""
        n = len(self.lat)
        if n == 0:
            return []
        span = float(self.width)
        if self.view_mode.active == 0:
            w = span / n
            return [(i * w, w) for i in range(n)]

        total = self.lat.total_length(self.lambda_w) or 1.0
        scale = span / total
        raw = [(zs * scale, L * scale) for _, zs, L in self.lat.layout(self.lambda_w)]
        # Zero-length elements (quads, modulations) still need to be
        # grabbable, so they get a floor width and everything after them
        # shifts along. The bar positions stay monotonic, which is what the
        # eye reads as "in order".
        out, shift = [], 0.0
        for x, w in raw:
            width = max(w, MIN_W)
            out.append((x + shift, width))
            shift += width - w
        squeeze = span / max(1e-9, out[-1][0] + out[-1][1]) if out else 1.0
        return [(x * squeeze, w * squeeze) for x, w in out]

    def _insertion_from_frac(self, frac):
        x = frac * self.width
        for i, (x0, w) in enumerate(self._slots_cache):
            if x < x0 + w / 2.0:
                return i
        return len(self.lat)

    # ── drawing ─────────────────────────────────────────────────────────

    def _strip_css(self):
        T = self.T
        return """
        <style>
          .pf-wrap { position:relative; width:100%%; }
          .pf-strip { display:block; width:100%%; height:%(h)dpx;
                      overflow:visible; touch-action:none; }
          .pf-bar { cursor:grab; transition:opacity .12s ease,
                    filter .12s ease, y .18s ease, height .18s ease; }
          .pf-bar:hover, .pf-bar.pf-hover { filter:brightness(1.12); }
          .pf-bar.pf-dragging { cursor:grabbing; opacity:.55;
                                filter:brightness(1.15); }
          .pf-sel { pointer-events:none; }
          .pf-drop { opacity:0; transition:opacity .1s ease;
                     pointer-events:none; }
          .pf-axis { font-family:%(font)s; font-size:%(micro)s;
                     fill:%(muted)s; }
          .pf-tick { stroke:%(grid)s; stroke-width:1; }
        </style>""" % dict(h=STRIP_H, font=pus.FONT, micro=pus.TYPE['micro'],
                           muted=T['muted'], grid=T['grid'])

    def _empty_strip(self, message):
        T = self.T
        return ('%s<div class="pf-wrap" style="height:%dpx; display:flex; '
                'align-items:center; justify-content:center; '
                'border:1px dashed %s; border-radius:%s;">%s</div>'
                % (self._strip_css(), STRIP_H, T['grid'], pus.RADIUS['md'],
                   pus.hint(T, message)))

    def redraw(self):
        if not self.enabled:
            return
        T = self.T
        slots = self._slots()
        self._slots_cache = slots
        if not slots:
            self.strip.text = self._empty_strip(
                'Empty lattice — add an element to start.')
            self._update_summary()
            self._build_fields()
            return

        W, gap = float(self.width), 1.2
        parts = [self._strip_css(),
                 '<div class="pf-wrap"><svg class="pf-strip" '
                 'viewBox="0 0 %g %d" preserveAspectRatio="none">' % (W, STRIP_H)]

        for i, (el, (x, w)) in enumerate(zip(self.lat.elements, slots)):
            h = BAR_H * HEIGHT.get(el.tag, 0.6)
            y = BAR_TOP + (BAR_H - h) / 2.0
            colour = T[COLOUR_ROLE.get(el.tag, 'drift')]
            bw = max(w - gap, w * 0.6)
            # Radius must not exceed half the bar, or a narrow bar turns
            # into a lozenge and stops reading as a block.
            r = min(2.5, bw / 2.0, h / 2.0)
            opacity = ' opacity="0.5"' if el.tag == 'DR' else ''
            parts.append(
                '<rect class="pf-bar" data-i="%d" x="%g" y="%g" width="%g" '
                'height="%g" rx="%g" fill="%s"%s><title>%s %d</title></rect>'
                % (i, x, y, bw, h, r, colour, opacity, el.label, i + 1))

        if self.selected is not None and 0 <= self.selected < len(slots):
            x, w = slots[self.selected]
            el = self.lat.elements[self.selected]
            h = BAR_H * HEIGHT.get(el.tag, 0.6)
            y = BAR_TOP + (BAR_H - h) / 2.0
            bw = max(w - gap, w * 0.6)
            parts.append(
                '<rect class="pf-sel" x="%g" y="%g" width="%g" height="%g" '
                'rx="%g" fill="none" stroke="%s" stroke-width="2"/>'
                % (x - 1.5, y - 3, bw + 3, h + 6,
                   min(4.0, (bw + 3) / 2.0), T['ink']))

        parts.append('<line class="pf-drop" x1="0" y1="%d" x2="0" y2="%d" '
                     'stroke="%s" stroke-width="2.5" stroke-linecap="round"/>'
                     % (BAR_TOP - 4, BAR_TOP + BAR_H + 4, T['ink']))
        parts.append(self._axis_svg(W))
        parts.append('</svg></div>')
        self.strip.text = ''.join(parts)

        self._update_summary()
        self._build_fields()

    def _axis_svg(self, W):
        """A light scale under the bars: z in metres, or element ordinals.

        The unit label sits at the right end, so ticks stop short of it —
        otherwise the last one prints straight through the label.
        """
        T = self.T
        y = BAR_TOP + BAR_H + 8
        label = 'element' if self.view_mode.active == 0 else 'z (m)'
        # Reserve room for the label: ~6.2px per character at 10px type.
        limit = W - (len(label) * 6.2 + 10)
        out = ['<line class="pf-tick" x1="0" y1="%g" x2="%g" y2="%g"/>'
               % (y, W, y)]

        def tick(x, text):
            if x > limit:
                return
            out.append('<text class="pf-axis" x="%g" y="%g" '
                       'text-anchor="%s">%s</text>'
                       % (x, y + 13, 'start' if x <= 1 else 'middle', text))

        n = len(self.lat)
        if self.view_mode.active == 0:
            step = max(1, int(round(n / 8.0)))
            for i in range(0, n + 1, step):
                tick((i / float(n)) * W if n else 0, str(i))
        else:
            total = self.lat.total_length(self.lambda_w) or 1.0
            nice = _nice_step(total)
            v = 0.0
            while v <= total + 1e-9:
                tick((v / total) * W, '%g' % round(v, 3))
                v += nice
        out.append('<text class="pf-axis" x="%g" y="%g" text-anchor="end" '
                   'style="font-style:italic;">%s</text>' % (W, y + 13, label))
        return ''.join(out)

    def _update_summary(self):
        T = self.T
        if not len(self.lat):
            self.summary.text = pus.hint(T, 'No elements.')
            return
        counts = {}
        for el in self.lat.elements:
            counts[el.tag] = counts.get(el.tag, 0) + 1
        chips = ' '.join(
            pus.chip(T, '%d %s%s' % (counts[t], SPEC[t][1].lower(),
                                     's' if counts[t] != 1 else ''),
                     COLOUR_ROLE.get(t, 'muted') if t != 'DR' else 'muted')
            for t in ORDER if t in counts)
        self.summary.text = (
            '<div style="display:flex; align-items:center; gap:%s; '
            'flex-wrap:wrap; padding-top:%s;">%s'
            '<span style="font-family:%s; font-size:%s; color:%s; '
            'font-variant-numeric:tabular-nums; margin-left:%s;">'
            '%s steps · %.2f m</span></div>'
            % (pus.SP[2], pus.SP[1], chips, pus.FONT, pus.TYPE['small'],
               T['muted'], pus.SP[2], '{:,}'.format(self.lat.total_steps()),
               self.lat.total_length(self.lambda_w)))

    # ── selection panel ─────────────────────────────────────────────────

    def _build_fields(self):
        T = self.T
        self.field_widgets = {}
        if self.selected is None or not (0 <= self.selected < len(self.lat)):
            self.sel_title.text = pus.hint(
                T, 'Click an element to edit it. Drag to move it along the sequence.')
            self.field_row.children = []
            return
        el = self.lat.elements[self.selected]
        self.sel_title.text = (
            '<div style="display:flex; align-items:center; gap:%s; '
            'font-family:%s; font-size:%s; color:%s; margin-top:%s;">%s'
            '<b style="color:%s;">%s</b>'
            '<span style="color:%s;">element %d of %d</span></div>'
            % (pus.SP[3], pus.FONT, pus.TYPE['body'], T['ink2'], pus.SP[2],
               pus.dot(T[COLOUR_ROLE.get(el.tag, 'drift')]), T['ink'],
               el.label, T['muted'], self.selected + 1, len(self.lat)))

        widgets = []
        for name in EDIT_FIELDS.get(el.tag, el.names):
            value = el.get(name, '')
            if name == 'type':
                opts = list(UND_TYPES)
                if value not in [v for v, _ in opts]:
                    opts.append((value, '%s (from file)' % value))
                w = Select(title='type', options=opts, value=value,
                           width=168, stylesheets=[self.icss])
            else:
                w = TextInput(title=name, value=str(value), width=92,
                              stylesheets=[self.icss])
                if FIELD_HELP.get(name):
                    w.description = FIELD_HELP[name]
            w.on_change('value', self._make_field_handler(name))
            self.field_widgets[name] = w
            widgets.append(w)
        self.field_row.children = widgets

    def _make_field_handler(self, name):
        def handler(attr, old, new):
            if self._suppress or self.selected is None:
                return
            self.lat.elements[self.selected].set(name, new)
            self._commit(redraw=name in ('periods', 'zbar', 'steps'))
        return handler

    # ── pointer events from the browser ─────────────────────────────────

    def _on_bus(self, attr, old, new):
        if not new or not self.enabled:
            return
        try:
            msg = json.loads(new)
        except ValueError:
            return
        action, i = msg.get('action'), msg.get('i')
        if not isinstance(i, int) or not (0 <= i < len(self.lat)):
            return
        if action == 'select':
            self.selected = i
            self.redraw()
        elif action == 'drop':
            target = self._insertion_from_frac(float(msg.get('f', 0.0)))
            # Dropping past its own slot shifts the target left by one,
            # because the element is removed before being reinserted.
            if target > i:
                target -= 1
            self.selected = self.lat.move(i, target)
            self._commit()

    # ── mutations ───────────────────────────────────────────────────────

    def _add(self):
        at = (self.selected + 1) if self.selected is not None else len(self.lat)
        self.selected = self.lat.insert(self.add_select.value, at)
        self._commit()

    def _duplicate(self):
        if self.selected is None:
            return
        self.selected = self.lat.duplicate(self.selected)
        self._commit()

    def _delete(self):
        if self.selected is None:
            return
        self.selected = self.lat.delete(self.selected)
        if not len(self.lat):
            self.selected = None
        self._commit()

    def _nudge(self, delta):
        if self.selected is None:
            return
        self.selected = self.lat.move(self.selected, self.selected + delta)
        self._commit()

    def _on_text_edit(self, attr, old, new):
        if self._suppress:
            return
        self.lat = Lattice.parse(new)
        if self.selected is not None and self.selected >= len(self.lat):
            self.selected = None
        self.redraw()
        self.on_change(new)

    def _commit(self, redraw=True):
        text = self.lat.text()
        self._suppress = True
        self.text.value = text
        self._suppress = False
        if redraw:
            self.redraw()
        else:
            self._update_summary()
        self.on_change(text)


def _nice_step(total):
    """A round tick interval giving roughly 6–10 ticks across `total`."""
    import math
    if total <= 0:
        return 1.0
    raw = total / 8.0
    mag = 10 ** math.floor(math.log10(raw))
    for mult in (1, 2, 2.5, 5, 10):
        if raw <= mult * mag:
            return mult * mag
    return 10 * mag
