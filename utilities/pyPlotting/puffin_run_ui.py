#!/usr/bin/env python3
"""Web UI for setting up and launching Puffin runs.

Run with:
    bokeh serve --show puffin_run_ui.py [--port 5685] [--args runs=/path/to/runs puffin=/path/to/puffin]

Three columns, left to right: pick an example deck, edit its parameters,
launch and watch. The flow is deliberately example-first — a Puffin deck has
enough interacting parameters that starting from a working one and changing
what you mean to change is the only sane way in, and it means the form never
has to invent defaults that the code would reject.

Editing is round-trip safe. The form promotes the parameters people tune
(see puffin_ui_fields.py); everything else in the deck keeps its value, its
comments and its layout, and stays editable in the raw tabs. A field is only
written back if the deck already carried it or you actually changed it, so
loading an example and pressing Go reproduces that example exactly rather
than a version silently padded with defaults.

Pressing Go copies the four input files into a fresh timestamped directory
and launches there, detached. Nothing runs in the browser and nothing is
overwritten: the run directory is the record of what was run.

Save only writes the same directory and stops, adding a `run.sh` that
carries the launch command — for a cluster, where the job belongs to a
scheduler rather than to this process. Both buttons go through the same
preparation, so a saved deck is byte-for-byte the one Go would have used.

Shares puffin_viz_theme with the field viewers, so the two read as one tool.
"""

import os
import sys
import glob
import time
import shutil
import datetime
import subprocess

from bokeh.plotting import curdoc
from bokeh.models import (Button, Div, Select, TextInput, NumericInput, Switch,
                          TextAreaInput, Tabs, TabPanel, InlineStyleSheet,
                          ScrollBox)
from bokeh.layouts import column, row

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import puffin_viz_theme as pvt
import puffin_ui_fields as puf
import puffin_ui_style as pus
import puffin_lattice_ui as plu
from puffin_deck import Deck, total_steps, format_value, parse_value
from puffin_lattice import Lattice
import puffin_runner as pr

T = pvt.tokens()
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))

POLL_MS = 700          # UI refresh while a run is live


# ── launch arguments ────────────────────────────────────────────────────────

def _args():
    out = {}
    for arg in sys.argv[1:]:
        if '=' in arg:
            key, val = arg.split('=', 1)
            out[key.strip()] = val.strip()
    return out


ARGS = _args()
RUNS_ROOT = os.path.abspath(os.path.expanduser(
    ARGS.get('runs', os.path.join(os.path.expanduser('~'), 'puffin-runs'))))
EXAMPLES_ROOT = os.path.abspath(os.path.expanduser(
    ARGS.get('examples', os.path.join(REPO, 'inputs'))))
PUFFIN_BIN = pr.find_puffin(ARGS.get('puffin'))


# ── examples ────────────────────────────────────────────────────────────────

def find_examples(root):
    """Directories holding a main deck, as {label: (dir, main filename)}.

    A main deck is the file carrying the &MDATA block; beam and seed files
    live beside it. Directories with more than one are listed once per deck.
    """
    found = {}
    for path in sorted(glob.glob(os.path.join(root, '**', '*.in'), recursive=True)):
        try:
            # Read the whole file: the &MDATA block sits below a documentation
            # header that runs to ~90 lines in most decks, so a head-only peek
            # misses it. Decks are a few KB.
            with open(path, 'r', errors='replace') as fh:
                text = fh.read()
        except OSError:
            continue
        if '&mdata' not in text.lower():
            continue
        directory, fname = os.path.split(path)
        label = os.path.relpath(directory, root)
        if len([p for p in found.values() if p[0] == directory]) or label in found:
            label = os.path.join(label, fname)
        found[label] = (directory, fname)
    return found


EXAMPLES = find_examples(EXAMPLES_ROOT)


# ── deck state ──────────────────────────────────────────────────────────────

class DeckSet:
    """The four files of one setup, and where they came from."""

    def __init__(self):
        self.source_dir = None
        self.main_name = None
        self.main = None
        self.beam = None
        self.seed = None
        self.beam_name = None
        self.seed_name = None
        self.lattice_name = None
        self.lattice_text = None

    def load(self, directory, main_name):
        self.source_dir = directory
        self.main_name = main_name
        self.main = Deck.load(os.path.join(directory, main_name))

        self.beam_name = self._referenced('beam_file')
        self.seed_name = self._referenced('seed_file')
        self.lattice_name = self._referenced('lattFile')

        self.beam = self._load_side(self.beam_name)
        self.seed = self._load_side(self.seed_name)
        self.lattice_text = None
        if self.lattice_name:
            path = os.path.join(directory, self.lattice_name)
            if os.path.isfile(path):
                with open(path, 'r', errors='replace') as fh:
                    self.lattice_text = fh.read()
        return self

    def _referenced(self, key):
        value = self.main.get(key, '')
        return value.strip() if isinstance(value, str) and value.strip() else None

    def _load_side(self, name):
        if not name:
            return None
        path = os.path.join(self.source_dir, name)
        return Deck.load(path) if os.path.isfile(path) else None

    def deck_for(self, which):
        return {'main': self.main, 'beam': self.beam, 'seed': self.seed}.get(which)

    def total_steps(self):
        """Steps the run will take, counted from the lattice as it stands now.

        Read from `lattice_text` rather than the file on disk: the editor
        works on the text, so a lattice you have just changed must be what
        the progress estimate is based on.
        """
        if self.lattice_text:
            return Lattice.parse(self.lattice_text).total_steps() or None
        if self.lattice_name:
            return None            # a lattice is named but could not be read
        return total_steps(self.main, None)

    def write_to(self, directory):
        """Write all four files into `directory`, copying any extra inputs."""
        os.makedirs(directory, exist_ok=True)
        self.main.save(os.path.join(directory, self.main_name))
        if self.beam and self.beam_name:
            self.beam.save(os.path.join(directory, self.beam_name))
        if self.seed and self.seed_name:
            self.seed.save(os.path.join(directory, self.seed_name))
        if self.lattice_name and self.lattice_text is not None:
            with open(os.path.join(directory, self.lattice_name), 'w') as fh:
                fh.write(self.lattice_text)
        # Anything else the deck leans on (dist files, h5 beams) travels too.
        written = {self.main_name, self.beam_name, self.seed_name, self.lattice_name}
        if self.source_dir:
            for name in os.listdir(self.source_dir):
                src = os.path.join(self.source_dir, name)
                if name in written or not os.path.isfile(src):
                    continue
                if name.endswith(('.in', '.latt', '.h5', '.pufin')):
                    shutil.copy(src, os.path.join(directory, name))
        return directory


DECKS = DeckSet()
RUNS = []           # newest first
ACTIVE = {'run': None}


# ── chrome ──────────────────────────────────────────────────────────────────

def _btn_css(kind='plain'):
    return pus.button_css(T, kind)


INPUT_CSS = pus.input_css(T)
TABS_CSS = pus.tabs_css(T)
FILL_CSS = pus.fill_css()


def section(title, note=''):
    """A titled group heading."""
    return Div(text=pus.heading(T, title, note),
               sizing_mode='stretch_width', stylesheets=[FILL_CSS])


def eyebrow(title):
    return Div(text=pus.eyebrow(T, title),
               sizing_mode='stretch_width', stylesheets=[FILL_CSS])


def hint(text):
    return Div(text=pus.hint(T, text),
               sizing_mode='stretch_width', stylesheets=[FILL_CSS])


# ── form ────────────────────────────────────────────────────────────────────

WIDGETS = {}        # (file, key) → widget
ORIGINAL = {}       # (file, key) → value text as loaded, or None if absent
SNAPSHOT = {}       # (file, key) → widget state right after populate
UNREPRESENTABLE = set()   # scalar fields holding a vector — raw tab only


def make_widget(field):
    kind, label = field['kind'], field['label']
    title = f"{label}  ({field['unit']})" if field['unit'] else label
    common = dict(title=title, stylesheets=[INPUT_CSS], width=190)
    if kind == 'bool':
        # Switch carries no title of its own, so it is labelled by the row.
        return Switch(active=False, width=42)
    if kind.startswith('choice:'):
        opts = puf.choices(kind)
        labels = puf.CHOICE_LABELS.get(field['key'], {})
        return Select(options=[(o, labels.get(o, o or '(none)')) for o in opts], **common)
    if kind == 'int':
        return NumericInput(mode='int', **common)
    if kind == 'float':
        return NumericInput(mode='float', **common)
    return TextInput(**common)


def build_form():
    """Group panels of widgets, one tab per group."""
    panels = []
    for group in puf.GROUPS:
        rows, pending, last_kind = [], [], None
        for field in group['fields']:
            key = (group['file'], field['key'])
            widget = make_widget(field)
            WIDGETS[key] = widget
            last_kind = field['kind']
            if field['kind'] == 'bool':
                # Switch takes neither a title nor a description, so the row
                # label carries both the name and the help (as a tooltip).
                tip = (' title="%s"' % field['help'].replace('"', '&quot;')
                       ) if field['help'] else ''
                lbl = Div(text=f"""<div{tip} style="font-family:{pvt.FONT};
                    font-size:12px; color:{T['ink']}; padding-top:3px;
                    {'border-bottom:1px dotted ' + T['axis'] + ';' if tip else ''}
                    display:inline-block;">{field['label']}</div>""", width=150)
                pending.append(row(widget, lbl, spacing=6))
                if len(pending) == 2:
                    rows.append(row(*pending, spacing=24)); pending = []
            else:
                if field['help']:
                    widget.description = field['help']
                pending.append(widget)
                if len(pending) == 2:
                    rows.append(row(*pending, spacing=16)); pending = []
        if pending:
            rows.append(row(*pending, spacing=24 if last_kind == 'bool' else 16))
        body = column(*rows, spacing=8, sizing_mode='stretch_width')
        panels.append(TabPanel(child=ScrollBox(child=body, sizing_mode='stretch_both'),
                               title=group['name']))
    return Tabs(tabs=panels, sizing_mode='stretch_both', stylesheets=[TABS_CSS])


def populate_form():
    """Load current deck values into the widgets, remembering what was there.

    Two things are recorded per field: the deck's own text (None when the key
    is absent) and the widget state that leaves this function. The second is
    what makes "absent and untouched" distinguishable from "deliberately
    set", because a widget has no empty state to fall back on — a Select
    lands on its first option whether or not the deck said anything.
    """
    ORIGINAL.clear()
    SNAPSHOT.clear()
    UNREPRESENTABLE.clear()
    for group, field in puf.all_fields():
        key = (group['file'], field['key'])
        widget = WIDGETS[key]
        deck = DECKS.deck_for(group['file'])
        raw = deck.raw(field['key']) if deck is not None else None
        ORIGINAL[key] = raw
        widget.disabled = deck is None
        if raw is None:
            _set_widget(widget, field, None)
        else:
            value = parse_value(raw) if field['kind'] != 'list' else raw
            # A multi-beam or multi-seed deck gives per-beam vectors to keys
            # this form treats as scalars — `emitx = 1.0, 1.0`, and even
            # `qFlatTop = .true., .true.` for a switch. One control cannot
            # stand for two beams, so the field steps aside and leaves it to
            # the raw tab, rather than showing the first entry and silently
            # collapsing the rest on write.
            if isinstance(value, list) and field['kind'] != 'list':
                UNREPRESENTABLE.add(key)
                widget.disabled = True
                _set_widget(widget, field, None)
            else:
                _set_widget(widget, field, value)
        SNAPSHOT[key] = _widget_state(widget, field)


def _widget_state(widget, field):
    """The widget's current value, for change detection."""
    return widget.active if field['kind'] == 'bool' else widget.value


def _same_value(new_text, old_text):
    """Do two namelist value texts mean the same thing?

    Compared as values, not strings: a deck writing `120` that comes back
    through a float box as `120.0` has not been edited, and reporting it as
    a change would bury the real ones in noise.
    """
    if new_text.strip() == old_text.strip():
        return True
    try:
        return parse_value(new_text) == parse_value(old_text)
    except Exception:
        return False


def _set_widget(widget, field, value):
    kind = field['kind']
    if kind == 'bool':
        widget.active = bool(value) if value is not None else False
    elif kind.startswith('choice:'):
        opts = puf.choices(kind)
        text = '' if value is None else str(value)
        widget.value = text if text in opts else (opts[0] if opts else '')
    elif kind in ('int', 'float'):
        widget.value = None if value in (None, '') else value
    else:
        widget.value = '' if value is None else str(value)


def _widget_text(widget, field):
    """Widget state → namelist value text, or None to leave the key alone."""
    kind = field['kind']
    if kind == 'bool':
        return '.true.' if widget.active else '.false.'
    if kind.startswith('choice:'):
        return format_value(widget.value) if widget.value != '' else "''"
    if kind in ('int', 'float'):
        if widget.value is None:
            return None
        return format_value(int(widget.value) if kind == 'int' else float(widget.value))
    text = (widget.value or '').strip()
    if not text:
        return None
    return text if kind == 'list' else format_value(text)


def collect_form():
    """Write changed / already-present fields back into the decks.

    A field the deck never carried, and that the user has not touched, stays
    absent: pressing Go on an unmodified example reproduces that example,
    not a version padded out with this form's idea of a default. Getting
    this wrong is not cosmetic — adding a key the deck's namelist group does
    not accept makes the file unreadable to Fortran.

    New keys go to the group's declared namelist block, never to whichever
    block happens to come first in the file.
    """
    changed = 0
    for group, field in puf.all_fields():
        key = (group['file'], field['key'])
        deck = DECKS.deck_for(group['file'])
        if deck is None:
            continue
        if key in UNREPRESENTABLE:
            continue                      # per-beam vector: raw tab owns it
        widget = WIDGETS[key]
        was = ORIGINAL.get(key)
        if was is None and _widget_state(widget, field) == SNAPSHOT.get(key):
            continue                      # absent from the deck, and untouched
        text = _widget_text(widget, field)
        if text is None:
            continue                      # cleared numeric: leave the deck alone
        if was is not None and _same_value(text, was):
            continue
        block = group.get('block')
        if was is None and block and block not in deck.blocks:
            # This deck uses a different namelist group for its beam (a
            # macroparticle or distribution deck carries &BDLIST, not
            # &BLIST), so the key has nowhere valid to go. Skipping is the
            # only safe move: appending it anywhere would break the read.
            continue
        deck.set(field['key'], text, block=block)
        changed += 1
    return changed


# ── raw editors ─────────────────────────────────────────────────────────────

RAW = {}
RAW_SNAPSHOT = {}   # text as last loaded, to tell a real edit from a redisplay


LATTICE = {'editor': None}


def _on_lattice_change(text):
    """The visual editor is the lattice's owner; keep the deck in step."""
    DECKS.lattice_text = text
    RAW_SNAPSHOT['lattice'] = text
    refresh_steps()


def build_raw():
    panels = []
    for which, label in (('main', 'Main deck'), ('beam', 'Beam'), ('seed', 'Seed')):
        area = TextAreaInput(value='', rows=26, sizing_mode='stretch_both',
                             stylesheets=[INPUT_CSS])
        RAW[which] = area
        panels.append(TabPanel(child=area, title=label))
    editor = plu.LatticeEditor(T, on_change=_on_lattice_change, width=940,
                               doc=curdoc())
    LATTICE['editor'] = editor
    RAW['lattice'] = editor.text          # the editor's own text box
    panels.append(TabPanel(child=editor.layout, title='Lattice'))
    return Tabs(tabs=panels, sizing_mode='stretch_both', stylesheets=[TABS_CSS])


def populate_raw():
    RAW_SNAPSHOT.clear()
    for which in ('main', 'beam', 'seed'):
        deck = DECKS.deck_for(which)
        RAW[which].value = deck.text() if deck is not None else ''
        RAW[which].disabled = deck is None
        RAW_SNAPSHOT[which] = RAW[which].value
    editor = LATTICE['editor']
    if editor is not None:
        has_lattice = DECKS.lattice_text is not None
        editor.set_enabled(has_lattice)
        if has_lattice:
            lw = DECKS.main.get('lambda_w') if DECKS.main is not None else None
            editor.set_text(DECKS.lattice_text,
                            lambda_w=lw if isinstance(lw, (int, float)) else None)
    RAW_SNAPSHOT['lattice'] = DECKS.lattice_text or ''


def apply_raw():
    """Adopt raw-tab edits. Returns True if any were found.

    A file counts as edited only when its text differs from what was last
    *displayed* there, not from the deck as it now stands — the form may
    have changed the deck since, and comparing against that would make every
    form edit look like a stale raw tab about to be undone.

    When the raw text really has been edited it wins for that file, and the
    form is refreshed from it. That direction is the safe one: raw is the
    escape hatch, so what you typed there is what runs.
    """
    edited = False
    for which in ('main', 'beam', 'seed'):
        deck = DECKS.deck_for(which)
        if deck is None:
            continue
        text = RAW[which].value
        if text.strip() and text != RAW_SNAPSHOT.get(which):
            setattr(DECKS, which, Deck.loads(text, path=deck.path))
            edited = True
    # The lattice editor writes straight through to DECKS.lattice_text as it
    # is used, so there is nothing to adopt here — only the three namelist
    # files have a raw tab that can drift from the model.

    if edited:
        populate_form()          # form edits are superseded by the raw text
        RAW_SNAPSHOT.update({k: RAW[k].value for k in RAW})
    return edited


# ── widgets: header, controls, run panel ────────────────────────────────────

example_select = Select(
    title='Example deck',
    options=sorted(EXAMPLES.keys()),
    value=(sorted(EXAMPLES.keys())[0] if EXAMPLES else ''),
    stylesheets=[INPUT_CSS], sizing_mode='stretch_width')

run_name = TextInput(title='Run name', value='', stylesheets=[INPUT_CSS], width=200)
# The name follows the selected deck until the user types one of their own —
# otherwise switching decks leaves the previous deck's name in place and the
# run directory ends up labelled after something it did not run.
NAME_EDITED = {'by_user': False}
ranks_input = NumericInput(title='MPI ranks', mode='int', value=2, low=1, high=256,
                           stylesheets=[INPUT_CSS], width=100)
ranks_input.description = (
    'Processes to run on. Beam loading and the parallel field decomposition '
    'both depend on this, so it is part of the physics setup, not just a '
    'speed knob — compare like with like when changing it.')

go_btn = Div()          # replaced below (needs handlers defined first)
status_div = Div(sizing_mode='stretch_width', stylesheets=[FILL_CSS])
progress_div = Div(sizing_mode='stretch_width', stylesheets=[FILL_CSS])
log_div = Div(sizing_mode='stretch_width', stylesheets=[FILL_CSS])
history_div = Div(sizing_mode='stretch_width', stylesheets=[FILL_CSS])
cmd_div = Div(sizing_mode='stretch_width', stylesheets=[FILL_CSS])


def progress_html(fraction, label, tone='field'):
    """A slim determinate bar, or a moving stripe when the total is unknown."""
    colour = T[tone]
    pct = 0 if fraction is None else int(round(100 * fraction))
    width = '100%' if fraction is None else '%d%%' % pct
    if fraction is None:
        fill = ('background-image:repeating-linear-gradient(115deg,'
                '%s 0 10px, %s33 10px 20px); '
                'animation:pfslide 1s linear infinite;' % (colour, colour))
    else:
        fill = 'background:%s;' % colour
    return """
    <style>@keyframes pfslide {{ to {{ background-position: 20px 0; }} }}</style>
    <div style="font-family:{font}; margin:{m} 0 0 0; width:100%; box-sizing:border-box;">
      <div style="height:6px; border-radius:999px; background:{grid};
                  overflow:hidden;">
        <div style="height:100%; width:{w}; {fill}
                    border-radius:999px;
                    transition:width .35s cubic-bezier(.4,0,.2,1);"></div>
      </div>
      <div style="font-size:{small}; color:{muted}; margin-top:{g};
                  font-variant-numeric:tabular-nums;">{label}</div>
    </div>""".format(font=pus.FONT, grid=T['grid'], w=width, fill=fill,
                     small=pus.TYPE['small'], muted=T['muted'],
                     m=pus.SP[3], g=pus.SP[3], label=label)


def status_html(text, tone='ink2'):
    return ('<div style="font-family:%s; font-size:%s; color:%s; '
            'line-height:1.55; padding:%s 0;">%s</div>'
            % (pus.FONT, pus.TYPE['small'], T[tone], pus.SP[1], text))


def log_html(text):
    body = (text or '').replace('&', '&amp;').replace('<', '&lt;')
    empty = ('<span style="color:%s;">no output yet</span>' % T['muted'])
    return ("""<div style="font-family:{mono}; font-size:{small};
        line-height:1.5; color:{ink2}; background:{surface};
        border:1px solid {grid}; border-radius:{r}; padding:{p} {px};
        height:180px; overflow:auto; white-space:pre-wrap;
        width:100%; box-sizing:border-box;
        box-shadow:inset 0 1px 2px rgba(0,0,0,.04);">{body}</div>"""
        .format(mono=pus.MONO, small=pus.TYPE['small'], ink2=T['ink2'],
                surface=T['surface'], grid=T['grid'], r=pus.RADIUS['md'],
                p=pus.SP[4], px=pus.SP[5], body=body or empty))


STATUS_TONE = {pr.STATUS_RUNNING: 'field', pr.STATUS_DONE: 'select',
               pr.STATUS_FAILED: 'beam', pr.STATUS_STOPPED: 'muted',
               pr.STATUS_READY: 'muted'}


# ── actions ─────────────────────────────────────────────────────────────────

def load_example(label):
    if label not in EXAMPLES:
        return
    directory, main_name = EXAMPLES[label]
    DECKS.load(directory, main_name)
    populate_form()
    populate_raw()
    if not NAME_EDITED['by_user']:
        auto = pr.slugify(os.path.basename(directory))
        if run_name.value != auto:
            run_name.value = auto
            NAME_EDITED['by_user'] = False      # our own write, not the user's
    refresh_command()
    refresh_steps()


def refresh_steps():
    """Redraw the deck summary line, including the current step count.

    Called on load and again whenever the lattice changes, so the step
    figure tracks edits instead of describing the deck as it arrived.
    """
    if DECKS.main is None:
        return
    steps = DECKS.total_steps()
    bits = [f'<b>{DECKS.main_name}</b>']
    bits.append(f'beam: {DECKS.beam_name}' if DECKS.beam_name else 'no beam file')
    if DECKS.seed_name:
        bits.append(f'seed: {DECKS.seed_name}')
    if DECKS.lattice_name:
        bits.append(f'lattice: {DECKS.lattice_name}')
    bits.append(f'{steps:,} steps' if steps else 'step count unknown')
    status_div.text = status_html(' &nbsp;·&nbsp; '.join(bits))


def refresh_command():
    if DECKS.main is None:
        cmd_div.text = ''
        return
    probe = pr.Run('<run dir>', DECKS.main_name, ranks=int(ranks_input.value or 1),
                   puffin_bin=PUFFIN_BIN)
    cmd_div.text = (
        '<div style="font-family:%s; font-size:%s; color:%s; '
        'padding:%s 0 0 0; width:100%%; box-sizing:border-box; '
        'overflow-x:auto; white-space:nowrap; '
        'scrollbar-width:thin;">%s</div>'
        % (pus.MONO, pus.TYPE['micro'], T['muted'], pus.SP[3],
           probe.command_text()))


def _prepare_directory():
    """Adopt every pending edit and write the deck to a fresh run directory.

    Returns (directory, changed_count), or (None, reason) if it could not.
    Shared by Go and Save so that a saved deck is byte-for-byte the one a
    run would have used — if the two paths could diverge, a deck submitted
    to a scheduler would not be the thing you tested here.
    """
    apply_raw()                  # raw edits first: they win over the form
    changed = collect_form()     # then the form, onto whatever raw left behind
    populate_raw()               # show the decks as they will actually be written

    stamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    directory = os.path.join(RUNS_ROOT, '%s_%s' % (stamp, pr.slugify(run_name.value)))
    try:
        DECKS.write_to(directory)
    except OSError as exc:
        return None, 'Could not write run directory: %s' % exc
    return directory, changed


def _changed_note(changed):
    if not changed:
        return ' · deck unmodified'
    return ' · %d parameter%s changed' % (changed, '' if changed == 1 else 's')


def on_save():
    """Write the deck without running it, for submission to a scheduler."""
    if DECKS.main is None:
        status_div.text = status_html('Load an example deck first.', 'beam')
        return

    directory, changed = _prepare_directory()
    if directory is None:
        status_div.text = status_html(changed, 'beam')
        return

    ranks = int(ranks_input.value or 1)
    script = _write_launch_script(directory, ranks)
    status_div.text = status_html(
        'Saved to <code>%s</code>%s<br>'
        'Inputs plus <code>%s</code> with the launch command — add your '
        'scheduler directives to that, or call it from a job script.'
        % (directory, _changed_note(changed), os.path.basename(script)))


def _write_launch_script(directory, ranks):
    """A minimal run.sh carrying the exact command this UI would have used.

    Deliberately no #SBATCH or #PBS header: those differ per site (queue,
    account, wall time, node layout), and a guessed one is worse than none —
    it looks authoritative and silently submits to the wrong place. This is
    the part that is the same everywhere; the site-specific part is yours.
    """
    probe = pr.Run(directory, DECKS.main_name, ranks=ranks, puffin_bin=PUFFIN_BIN)
    path = os.path.join(directory, 'run.sh')
    threads = probe.environment_overrides().get('OMP_NUM_THREADS', '1')

    # The detected paths go in as defaults, not as literals. They are right
    # when the UI runs on the same machine as the job, and a `module load`
    # or a different compute-node image is exactly when they stop being
    # right — at which point overriding a variable beats editing the script.
    binary = PUFFIN_BIN or 'puffin'
    launcher = pr.find_mpiexec() or 'mpiexec'
    if ranks > 1:
        command = '"$MPIEXEC" -n "$RANKS" "$PUFFIN" %s' % DECKS.main_name
    else:
        command = '"$PUFFIN" %s' % DECKS.main_name

    with open(path, 'w') as fh:
        fh.write(
            '#!/bin/sh\n'
            '# Puffin run, written by puffin_run_ui.py on %s\n'
            '#\n'
            '# Add scheduler directives above, or call this from a job script.\n'
            '# The paths below are what this machine had when the deck was\n'
            '# saved; override them in the environment if the job runs\n'
            '# somewhere else.\n'
            '\n'
            'cd "$(dirname "$0")" || exit 1\n'
            '\n'
            ': "${PUFFIN:=%s}"\n'
            ': "${MPIEXEC:=%s}"\n'
            ': "${RANKS:=%d}"\n'
            '\n'
            '# Threading is a net loss with the default libgomp build, and\n'
            '# leaving this unset lets OpenMP claim cores on top of the MPI\n'
            '# ranks, which oversubscribes badly.\n'
            ': "${OMP_NUM_THREADS:=%s}"\n'
            'export OMP_NUM_THREADS\n'
            '\n'
            '%s\n'
            % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
               binary, launcher, ranks, threads, command))
    os.chmod(path, 0o755)
    return path


def on_go():
    if DECKS.main is None:
        status_div.text = status_html('Load an example deck first.', 'beam')
        return
    if ACTIVE['run'] is not None and ACTIVE['run'].status == pr.STATUS_RUNNING:
        status_div.text = status_html('A run is already in progress.', 'beam')
        return
    if not PUFFIN_BIN:
        status_div.text = status_html(
            'No puffin binary found — build it, or pass '
            '<code>--args puffin=/path/to/puffin</code>.', 'beam')
        return

    directory, changed = _prepare_directory()
    if directory is None:
        status_div.text = status_html(changed, 'beam')
        return

    steps = DECKS.total_steps()
    cadence = DECKS.main.get('iWriteIntNthSteps')
    run = pr.Run(directory, DECKS.main_name, ranks=int(ranks_input.value or 1),
                 puffin_bin=PUFFIN_BIN, expected_steps=steps,
                 write_int_every=cadence if isinstance(cadence, int) else None)
    try:
        run.start()
    except (RuntimeError, OSError) as exc:
        status_div.text = status_html('Could not start: %s' % exc, 'beam')
        return

    ACTIVE['run'] = run
    RUNS.insert(0, run)
    status_div.text = status_html(
        'Running in <code>%s</code>%s' % (directory, _changed_note(changed)))
    go_button.disabled = True
    stop_button.disabled = False
    view_button.disabled = True
    tick()


def on_stop():
    run = ACTIVE['run']
    if run is not None:
        run.stop()
        status_div.text = status_html('Stopped.', 'beam')
    tick()


def on_view():
    """Open the matching field viewer against the finished run."""
    run = ACTIVE['run']
    if run is None:
        return
    if not run.field_dumps():
        status_div.text = status_html(
            'No field dumps to view — this deck wrote none '
            '(iWriteNthSteps may exceed the run length).', 'beam')
        return
    viewer = 'viewField1D_bokeh.py' if run.is_1d(DECKS.main) else 'viewField3D_bokeh.py'
    base = os.path.splitext(DECKS.main_name)[0]
    port = _free_port()
    argv = [_bokeh_bin(), 'serve', '--show', os.path.join(HERE, viewer),
            '--port', str(port),
            '--allow-websocket-origin=localhost:%d' % port,
            '--args', run.dir, base]
    if DECKS.lattice_name:
        latt = os.path.join(run.dir, DECKS.lattice_name)
        if os.path.isfile(latt):
            argv.append(latt)
    try:
        subprocess.Popen(argv, cwd=REPO, start_new_session=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as exc:
        status_div.text = status_html('Could not launch viewer: %s' % exc, 'beam')
        return
    status_div.text = status_html(
        '%s starting on <a href="http://localhost:%d" target="_blank" '
        'style="color:%s">localhost:%d</a> — give it a moment.'
        % (viewer, port, T['field'], port))


def _bokeh_bin():
    """The bokeh launcher beside this interpreter, so the viewer runs in the
    same environment as the UI rather than whatever is first on PATH.

    Windows venvs put it in Scripts/ with an .exe suffix; POSIX ones in bin/.
    """
    base = os.path.dirname(sys.executable)
    for cand in (os.path.join(base, 'bokeh'),
                 os.path.join(base, 'bokeh.exe'),
                 os.path.join(base, 'Scripts', 'bokeh.exe')):
        if os.path.isfile(cand):
            return cand
    return shutil.which('bokeh') or 'bokeh'


def _free_port():
    import socket
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]


# ── live refresh ────────────────────────────────────────────────────────────

def fmt_dt(seconds):
    seconds = int(seconds)
    return '%d:%02d' % (seconds // 60, seconds % 60) if seconds >= 60 else '%ds' % seconds


def tick():
    run = ACTIVE['run']
    if run is None:
        progress_div.text = progress_html(0.0, 'idle')
        log_div.text = log_html('')
        return

    state = run.poll()
    fraction = run.progress()
    dumps = run.dump_count()
    module = run.current_module()

    parts = []
    if fraction is not None:
        parts.append('%d%%' % int(round(100 * fraction)))
    parts.append('%d dump%s' % (dumps, '' if dumps == 1 else 's'))
    if run.expected_dumps:
        parts[-1] += ' of ~%d' % run.expected_dumps
    if module:
        parts.append('module %d' % module)
    parts.append(fmt_dt(run.elapsed()))
    label = '%s · %s' % (state, ' · '.join(parts))

    progress_div.text = progress_html(
        fraction if state == pr.STATUS_RUNNING or fraction else
        (1.0 if state == pr.STATUS_DONE else fraction), label)
    log_div.text = log_html(run.log_tail(4000))
    history_div.text = history_html()

    if state != pr.STATUS_RUNNING:
        go_button.disabled = False
        stop_button.disabled = True
        view_button.disabled = not run.field_dumps()
        if state == pr.STATUS_DONE:
            status_div.text = status_html(
                'Finished in %s · <code>%s</code>' % (fmt_dt(run.elapsed()), run.dir))
        elif state == pr.STATUS_FAILED:
            status_div.text = status_html(
                'Failed: %s · <code>%s</code>' % (run.error or 'see run.log', run.dir), 'beam')


def history_html():
    if not RUNS:
        return pus.hint(T, 'Runs from this session will appear here.')
    rows = []
    for run in RUNS[:12]:
        tone = STATUS_TONE.get(run.status, 'muted')
        pulse = (' animation:pfpulse 1.4s ease-in-out infinite;'
                 if run.status == pr.STATUS_RUNNING else '')
        rows.append(
            '<div style="display:flex; align-items:center; gap:%s; '
            'padding:%s 0; border-bottom:1px solid %s;">'
            '<span style="width:7px;height:7px;border-radius:50%%;'
            'background:%s;flex:none;%s"></span>'
            '<span style="font-family:%s; font-size:%s; color:%s; flex:1; '
            'overflow:hidden; text-overflow:ellipsis; white-space:nowrap;">'
            '%s</span>'
            '<span style="font-family:%s; font-size:%s; color:%s; '
            'font-variant-numeric:tabular-nums;">%s</span></div>'
            % (pus.SP[3], pus.SP[2], T['grid'], T[tone], pulse,
               pus.FONT, pus.TYPE['small'], T['ink2'], run.name,
               pus.FONT, pus.TYPE['micro'], T['muted'], fmt_dt(run.elapsed())))
    # The keyframes block is kept out of the %-formatted string: `%{` in
    # `50%{...}` reads as a format spec and raises.
    css = '<style>@keyframes pfpulse{50%{opacity:.35}}</style>'
    return css + ('<div style="margin-top:%s;">%s</div>'
                  % (pus.SP[2], ''.join(rows)))


# ── assemble ────────────────────────────────────────────────────────────────

go_button = Button(label='Go', button_type='primary', width=90,
                   stylesheets=[_btn_css('go')])
stop_button = Button(label='Stop', width=80, disabled=True,
                     stylesheets=[_btn_css('stop')])
# Button carries no `description` in Bokeh 3.4, so the explanation lives in
# the hint under the row rather than in a tooltip.
save_button = Button(label='Save only', width=98,
                     stylesheets=[_btn_css('plain')])
view_button = Button(label='Open viewer', width=120, disabled=True,
                     stylesheets=[_btn_css('plain')])

go_button.on_click(on_go)
stop_button.on_click(on_stop)
save_button.on_click(on_save)
view_button.on_click(on_view)
example_select.on_change('value', lambda a, o, n: load_example(n))
ranks_input.on_change('value', lambda a, o, n: refresh_command())


def _on_name_change(attr, old, new):
    # load_example resets the flag straight after its own write, so only a
    # change this callback does not recognise as ours counts as the user's.
    if new.strip() and new != pr.slugify(os.path.basename(
            EXAMPLES.get(example_select.value, ('', ''))[0] or '')):
        NAME_EDITED['by_user'] = True


run_name.on_change('value', _on_name_change)

form_tabs = build_form()
raw_tabs = build_raw()

header = Div(text="""
<div style="font-family:%s; padding:0 0 %s 0;">
  <div style="display:flex; align-items:center; gap:%s;">
    <svg width="26" height="26" viewBox="0 0 24 24" fill="none"
         stroke="%s" stroke-width="2" stroke-linecap="round"
         stroke-linejoin="round" aria-hidden="true">
      <path d="M2 13h3.5l2.5-7 3.5 13 2.5-9 2 3H22"/>
    </svg>
    <div>
      <div style="font-size:%s; letter-spacing:.11em; text-transform:uppercase;
                  color:%s; font-weight:600; line-height:1;">Puffin</div>
      <div style="font-size:%s; font-weight:600; color:%s;
                  line-height:1.15; margin-top:2px;">Run setup</div>
    </div>
  </div>
</div>""" % (pus.FONT, pus.SP[4], pus.SP[4], T['field'], pus.TYPE['micro'],
             T['muted'], pus.TYPE['title'], T['ink']),
    sizing_mode='stretch_width')


def card(*children, **kw):
    """A raised panel. Bokeh has no card model, so this is a styled column."""
    return column(*children, sizing_mode=kw.pop('sizing_mode', 'stretch_width'),
                  spacing=kw.pop('spacing', 6),
                  styles={'background': T['surface'],
                          'border': '1px solid ' + T['grid'],
                          'border-radius': pus.RADIUS['lg'],
                          'padding': pus.SP[6],
                          'box-shadow': pus.shadow(T, 1)}, **kw)


left = column(
    header,
    card(eyebrow('Deck'),
         example_select,
         hint('Loads the deck, its beam and seed files, and any lattice.'),
         status_div),
    card(eyebrow('Run'),
         row(run_name, ranks_input, spacing=10),
         row(go_button, stop_button, save_button, spacing=6),
         row(view_button, spacing=6),
         hint('Save only writes the deck and a run.sh without launching — '
              'for submitting to a scheduler yourself.'),
         cmd_div),
    card(eyebrow('Progress'),
         progress_div),
    card(eyebrow('Output'),
         log_div),
    card(eyebrow('History'),
         history_div),
    width=440, spacing=12,
)

right = column(
    card(section('Parameters', 'the ones worth tuning'),
         form_tabs,
         sizing_mode='stretch_both'),
    card(section('Files', 'everything, verbatim'),
         hint('Adopted when you press Go. Editing a file here overrides the '
              'form above for that file; otherwise the form wins.'),
         raw_tabs,
         sizing_mode='stretch_both'),
    sizing_mode='stretch_both', spacing=12,
)

layout = row(left, right, sizing_mode='stretch_both', spacing=20,
             stylesheets=[pus.global_css(T)],
             styles={'padding': pus.SP[7] + ' ' + pus.SP[7] + ' ' +
                                pus.SP[8] + ' ' + pus.SP[7]})

doc = curdoc()
doc.title = 'Puffin — run setup'
pus.apply(doc, T)
doc.add_root(layout)

if EXAMPLES:
    load_example(example_select.value)
else:
    status_div.text = status_html(
        'No example decks found under <code>%s</code>. '
        'Pass <code>--args examples=/path/to/inputs</code>.' % EXAMPLES_ROOT, 'beam')
if not PUFFIN_BIN:
    cmd_div.text = status_html(
        'No puffin binary found — build it, or pass '
        '<code>--args puffin=/path/to/puffin</code>.', 'beam')

progress_div.text = progress_html(0.0, 'idle')
log_div.text = log_html('')
doc.add_periodic_callback(tick, POLL_MS)
