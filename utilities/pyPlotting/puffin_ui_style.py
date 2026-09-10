#!/usr/bin/env python3
"""Design system for the Puffin UIs: scale, chrome, and widget stylesheets.

Colour still comes from puffin_viz_theme — that module owns the palette and
its colour-vision validation, and nothing here second-guesses it. What this
adds is everything else a page needs to look deliberate rather than
assembled: a type scale, a spacing scale, radii, elevation, and the CSS that
carries them into Bokeh's widgets.

That last part is the awkward one. Bokeh renders each widget into its own
shadow root and sets its own font inside, so page-level CSS does not inherit
in — styling a widget means handing it a stylesheet of its own. The helpers
here exist so that happens once, consistently, instead of every module
inventing its own padding.

Type and space are scales, not free numbers. Picking from a small set is
what makes spacing look intentional; the moment sizes are chosen ad hoc per
element, the rhythm goes and no amount of individual tuning brings it back.
"""

import puffin_viz_theme as pvt

# ── scales ──────────────────────────────────────────────────────────────────

# Type. Steps are roughly 1.2x, rounded to whole pixels so text stays crisp.
TYPE = {
    'micro': '10px',   # axis ticks, unit suffixes
    'small': '11px',   # secondary text, help, labels
    'body': '12px',    # widget text, default
    'medium': '13px',  # emphasised body
    'large': '15px',   # card titles
    'title': '19px',   # section heading
    'hero': '26px',    # page title, the z readout
}

# Space. One doubling scale; anything not on it looks like a mistake.
SP = {0: '0', 1: '2px', 2: '4px', 3: '6px', 4: '8px',
      5: '12px', 6: '16px', 7: '24px', 8: '32px'}

RADIUS = {'sm': '6px', 'md': '10px', 'lg': '14px', 'pill': '999px'}

# Elevation. Two levels only: resting and raised. More than that and the
# page starts to read as a pile rather than a layout.
def shadow(t, level=1):
    if t['mode'] == 'dark':
        return ('0 1px 2px rgba(0,0,0,.5)' if level == 1
                else '0 4px 14px rgba(0,0,0,.55)')
    return ('0 1px 2px rgba(16,15,12,.06)' if level == 1
            else '0 4px 14px rgba(16,15,12,.10)')


# A webfont with a full system fallback. Inter if it loads, the platform UI
# face if not — which is what happens on an air-gapped cluster, and the
# reason the stack has to stand on its own rather than assume the network.
FONT_URL = 'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&display=swap'
FONT = pvt.FONT          # one stack, owned by the theme module
MONO = ("ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, "
        "'Liberation Mono', monospace")


# ── page shell ──────────────────────────────────────────────────────────────

def page_template(t):
    """The document shell: fonts, page ground, and the tokens as CSS vars.

    Bokeh prepends `{% extends base %}` itself, so this must not carry one.
    """
    return """{%% block preamble %%}
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="%(font_url)s" rel="stylesheet">
<style>
  :root {
    --pf-surface: %(surface)s;  --pf-page: %(page)s;
    --pf-ink: %(ink)s;          --pf-ink2: %(ink2)s;
    --pf-muted: %(muted)s;      --pf-grid: %(grid)s;
    --pf-axis: %(axis)s;        --pf-field: %(field)s;
    --pf-beam: %(beam)s;        --pf-select: %(select)s;
    --pf-radius: %(radius)s;    --pf-shadow: %(shadow)s;
  }
  html, body {
    background: %(page)s; margin: 0; padding: 0 0 %(pad)s 0;
    font-family: %(font)s;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
  }
  ::selection { background: %(field)s33; }
  /* Scrollbars, so long panes do not punch a bright hole in a dark page. */
  * { scrollbar-width: thin; scrollbar-color: %(axis)s transparent; }
  *::-webkit-scrollbar { width: 9px; height: 9px; }
  *::-webkit-scrollbar-thumb {
    background: %(axis)s; border-radius: 999px;
    border: 2px solid transparent; background-clip: content-box; }
  *::-webkit-scrollbar-track { background: transparent; }
</style>
{%% endblock %%}
""" % dict(font_url=FONT_URL, font=FONT, surface=t['surface'], page=t['page'],
           ink=t['ink'], ink2=t['ink2'], muted=t['muted'], grid=t['grid'],
           axis=t['axis'], field=t['field'], beam=t['beam'], select=t['select'],
           radius=RADIUS['md'], shadow=shadow(t, 1), pad=SP[8])


def apply(doc, t):
    doc.template = page_template(t)


# ── widget stylesheets ──────────────────────────────────────────────────────

def input_css(t):
    """Text, numeric, select and textarea chrome.

    Bokeh's own font does not reach past the shadow boundary, so the family
    is set explicitly here — without it the widgets sit in Helvetica while
    the page around them is in Inter.
    """
    from bokeh.models import InlineStyleSheet
    return InlineStyleSheet(css="""
      .bk-input {
        background: %(surface)s; color: %(ink)s;
        border: 1px solid %(axis)s; border-radius: %(r)s;
        font-family: %(font)s; font-size: %(body)s;
        padding: %(pad)s %(padx)s; line-height: 1.35;
        transition: border-color .14s ease, box-shadow .14s ease;
        box-shadow: none;
      }
      .bk-input:hover:not(:focus):not(:disabled) { border-color: %(muted)s; }
      .bk-input:focus {
        border-color: %(field)s; outline: none;
        box-shadow: 0 0 0 3px %(field)s26;
      }
      .bk-input:disabled { opacity: .5; cursor: not-allowed; }
      .bk-input-group > label {
        color: %(ink2)s; font-family: %(font)s;
        font-size: %(small)s; font-weight: 500;
        margin-bottom: %(gap)s; display: block;
        letter-spacing: .01em;
      }
      select.bk-input { cursor: pointer; }
      textarea.bk-input { font-family: %(mono)s; font-size: %(small)s;
                          line-height: 1.5; }
    """ % dict(surface=t['surface'], ink=t['ink'], ink2=t['ink2'],
               axis=t['axis'], muted=t['muted'], field=t['field'],
               font=FONT, mono=MONO, body=TYPE['body'], small=TYPE['small'],
               r=RADIUS['sm'], pad=SP[2], padx=SP[3], gap=SP[1]))


def button_css(t, kind='plain'):
    """Button chrome. `kind` picks the role: go, stop, plain, ghost."""
    from bokeh.models import InlineStyleSheet
    roles = {
        'go':    (t['field'],   '#ffffff',  t['field']),
        'stop':  (t['beam'],    '#ffffff',  t['beam']),
        'plain': (t['surface'], t['ink2'],  t['axis']),
        'ghost': ('transparent', t['muted'], 'transparent'),
    }
    bg, fg, border = roles.get(kind, roles['plain'])
    lift = shadow(t, 1) if kind in ('go', 'stop') else 'none'
    return InlineStyleSheet(css="""
      .bk-btn {
        background: %(bg)s; color: %(fg)s;
        border: 1px solid %(border)s; border-radius: %(r)s;
        font-family: %(font)s; font-size: %(body)s; font-weight: 500;
        padding: %(pady)s %(padx)s; letter-spacing: .005em;
        box-shadow: %(lift)s;
        transition: filter .14s ease, transform .1s ease,
                    box-shadow .14s ease, border-color .14s ease;
      }
      .bk-btn:hover:not(:disabled) {
        filter: brightness(%(hover)s);
        border-color: %(hover_border)s;
      }
      .bk-btn:active:not(:disabled) { transform: translateY(1px); }
      .bk-btn:focus-visible {
        outline: none; box-shadow: 0 0 0 3px %(field)s33;
      }
      .bk-btn:disabled { opacity: .4; box-shadow: none; cursor: not-allowed; }
    """ % dict(bg=bg, fg=fg, border=border, font=FONT, body=TYPE['body'],
               r=RADIUS['sm'], pady=SP[3], padx=SP[5], lift=lift,
               field=t['field'],
               hover='1.07' if kind in ('go', 'stop') else '1.0',
               hover_border=t['muted'] if kind in ('plain', 'ghost') else border))


def tabs_css(t):
    """Tab strip: a quiet underline rather than boxed folders."""
    from bokeh.models import InlineStyleSheet
    return InlineStyleSheet(css="""
      .bk-header { border-bottom: 1px solid %(grid)s; gap: %(gap)s; }
      .bk-tab {
        font-family: %(font)s; font-size: %(body)s; font-weight: 500;
        color: %(muted)s; background: transparent; border: none;
        padding: %(pady)s %(padx)s; border-radius: 0;
        border-bottom: 2px solid transparent; margin-bottom: -1px;
        transition: color .14s ease, border-color .14s ease;
      }
      .bk-tab:hover { color: %(ink2)s; }
      .bk-tab.bk-active {
        color: %(field)s; border-bottom-color: %(field)s;
        background: transparent;
      }
    """ % dict(grid=t['grid'], muted=t['muted'], ink2=t['ink2'],
               field=t['field'], font=FONT, body=TYPE['body'],
               pady=SP[3], padx=SP[4], gap=SP[2]))


def radio_css(t):
    """Segmented control for the small either/or choices."""
    from bokeh.models import InlineStyleSheet
    return InlineStyleSheet(css="""
      .bk-btn-group {
        background: %(grid)s; border-radius: %(r)s; padding: 2px; gap: 2px;
      }
      .bk-btn-group > .bk-btn {
        font-family: %(font)s; font-size: %(small)s; font-weight: 500;
        color: %(ink2)s; background: transparent; border: none;
        border-radius: %(ri)s; padding: %(pady)s %(padx)s;
        transition: background .14s ease, color .14s ease;
      }
      .bk-btn-group > .bk-btn:hover:not(.bk-active) { color: %(ink)s; }
      .bk-btn-group > .bk-btn.bk-active {
        background: %(surface)s; color: %(ink)s; box-shadow: %(sh)s;
      }
    """ % dict(grid=t['grid'], surface=t['surface'], ink=t['ink'],
               ink2=t['ink2'], font=FONT, small=TYPE['small'],
               r=RADIUS['sm'], ri='4px', pady=SP[2], padx=SP[4],
               sh=shadow(t, 1)))


def fill_css():
    """Make a Div's content fill the widget instead of shrink-wrapping.

    Bokeh wraps Div content in `.bk-clearfix`, which is `inline-block`. Any
    `width:100%` inside therefore resolves against a shrink-wrapped box
    rather than the widget, so a log pane or a progress bar collapses to the
    width of its own text. The host element is sized correctly; only this
    wrapper is in the way.
    """
    from bokeh.models import InlineStyleSheet
    return InlineStyleSheet(css="""
      .bk-clearfix { display: block; width: 100%; box-sizing: border-box; }
    """)


def global_css(t):
    """Rules that must reach every widget, whatever its own sheet says."""
    from bokeh.models import GlobalInlineStyleSheet
    return GlobalInlineStyleSheet(css="""
      :host, .bk-Row, .bk-Column { font-family: %(font)s; }
      .bk-tooltip-content, .bk-Tooltip {
        font-family: %(font)s !important; font-size: %(small)s !important;
        background: %(surface)s !important; color: %(ink2)s !important;
        border: 1px solid %(grid)s !important;
        border-radius: %(r)s !important; box-shadow: %(sh)s !important;
        padding: %(pad)s %(padx)s !important; max-width: 280px;
        line-height: 1.45 !important;
      }
    """ % dict(font=FONT, small=TYPE['small'], surface=t['surface'],
               ink2=t['ink2'], grid=t['grid'], r=RADIUS['sm'],
               sh=shadow(t, 2), pad=SP[3], padx=SP[4]))


# ── HTML fragments ──────────────────────────────────────────────────────────

def card_open(t, pad=None):
    return ('<div style="background:%s; border:1px solid %s; '
            'border-radius:%s; padding:%s; box-shadow:%s;">'
            % (t['surface'], t['grid'], RADIUS['lg'], pad or SP[6],
               shadow(t, 1)))


def eyebrow(t, text):
    """A small uppercase label above a group."""
    return ('<div style="font-family:%s; font-size:%s; font-weight:600; '
            'letter-spacing:.09em; text-transform:uppercase; color:%s; '
            'margin:0 0 %s 0;">%s</div>'
            % (FONT, TYPE['micro'], t['muted'], SP[3], text))


def heading(t, text, note=''):
    extra = ('<span style="font-weight:400; color:%s; text-transform:none; '
             'letter-spacing:0; font-size:%s;"> %s</span>'
             % (t['muted'], TYPE['small'], note)) if note else ''
    return ('<div style="font-family:%s; font-size:%s; font-weight:600; '
            'color:%s; margin:0 0 %s 0;">%s%s</div>'
            % (FONT, TYPE['large'], t['ink'], SP[3], text, extra))


def hint(t, text):
    return ('<div style="font-family:%s; font-size:%s; color:%s; '
            'line-height:1.5; margin:%s 0 0 0;">%s</div>'
            % (FONT, TYPE['small'], t['muted'], SP[1], text))


def chip(t, text, role='muted'):
    """A small status pill. `role` maps to a theme colour."""
    colour = {'muted': t['muted'], 'field': t['field'], 'beam': t['beam'],
              'select': t['select'], 'ink2': t['ink2']}.get(role, t['muted'])
    return ('<span style="display:inline-flex; align-items:center; gap:%s; '
            'font-family:%s; font-size:%s; font-weight:500; color:%s; '
            'background:%s1a; border:1px solid %s33; border-radius:%s; '
            'padding:2px %s; letter-spacing:.01em; white-space:nowrap;">'
            '%s</span>'
            % (SP[2], FONT, TYPE['micro'], colour, colour, colour,
               RADIUS['pill'], SP[3], text))


def dot(colour, size=7):
    return ('<span style="display:inline-block; width:%dpx; height:%dpx; '
            'border-radius:50%%; background:%s; flex:none;"></span>'
            % (size, size, colour))


def mono(t, text, size=None):
    return ('<span style="font-family:%s; font-size:%s; color:%s;">%s</span>'
            % (MONO, size or TYPE['small'], t['ink2'], text))
