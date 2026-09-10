#!/usr/bin/env python3
"""Which deck parameters get a first-class control, and how they present.

The full namelist runs to ~90 keys across three files. Putting all of them
on screen as a flat list would be complete but unreadable, so this is a
curated subset: the parameters that get tuned between runs, grouped the way
the physics groups them, with the units and help text that the deck headers
document. Everything omitted here stays reachable in the raw editor tab, and
round-trips untouched — curation decides what is *promoted*, never what is
allowed.

Each field is (key, label, kind, unit, help). `kind` drives the widget and
the parse on the way back:

    float int      numeric entry
    bool           switch, written as .true./.false.
    text           free text, quoted on write if it is not numeric
    choice:a|b|c   dropdown
    list           comma-separated vector, written back verbatim

`file` on a group says which deck the keys belong to: 'main', 'beam' or
'seed'. `block` names the namelist group a *new* key must be appended to —
it matters because the beam and seed files open with a small header block
(`&NBLIST`, `&NSLIST`) that accepts only a couple of names, so a key added
to the first block of the file rather than the right one makes the deck
unreadable to Fortran.
"""

Field = lambda key, label, kind, unit='', help='': dict(
    key=key, label=label, kind=kind, unit=unit, help=help)


GROUPS = [
    dict(name='Undulator & FEL', file='main', block='mdata', fields=[
        Field('srho', 'ρ (Pierce)', 'float', '',
              'FEL parameter — strength of the interaction, and the efficiency scale'),
        Field('saw', 'a_w', 'float', '', 'Peak undulator parameter'),
        Field('sgamma_r', 'γ_r', 'float', '', 'Resonant / reference beam energy'),
        Field('lambda_w', 'λ_w', 'float', 'm', 'Undulator period'),
        Field('zundType', 'Undulator type', 'choice:|planepole|curved|helical',
              '', 'Blank means 1D field with no off-axis variation of a_w'),
        Field('taper', 'Taper', 'float', 'd(a_w)/dz',
              'Gradient of the undulator taper'),
        Field('sux', 'u_x', 'float', '', 'Polarisation: 1 with u_y=1 is helical'),
        Field('suy', 'u_y', 'float', '', 'Polarisation: u_y=0 with u_x=1 is planar'),
    ]),

    dict(name='Field mesh', file='main', block='mdata', fields=[
        Field('meshType', 'Mesh', 'choice:0|1', '',
              '0 = temporal (full pulse), 1 = periodic (single slice)'),
        Field('iNumNodesX', 'Nodes x', 'int', '', 'Field sampling nodes in x'),
        Field('iNumNodesY', 'Nodes y', 'int', '', 'Field sampling nodes in y'),
        Field('nodesPerLambdar', 'Nodes per λ_r', 'int', '',
              'Longitudinal sampling — nodes per resonant wavelength'),
        Field('sFModelLengthX', 'Model length x', 'float', 'm', ''),
        Field('sFModelLengthY', 'Model length y', 'float', 'm', ''),
        Field('sFModelLengthZ2', 'Model length z2', 'float', '',
              'Length of the field model in scaled z2'),
        Field('sPerWaves', 'Periodic waves', 'float', '',
              'Radiation periods held in the mesh, periodic mode'),
        Field('iRedNodesX', 'Inner nodes x', 'int', '',
              'Central region in x the electrons are held within'),
        Field('iRedNodesY', 'Inner nodes y', 'int', ''),
        Field('sFiltFrac', 'Filter cutoff', 'float', '× ω_r',
              'High-pass cutoff as a fraction of the resonant frequency'),
        Field('sDiffFrac', 'Diffraction step', 'float', '× λ_w',
              'Diffraction step size as a fraction of the undulator period'),
    ]),

    dict(name='Integration & output', file='main', block='mdata', fields=[
        Field('lattFile', 'Lattice file', 'text', '',
              'Optional. When set, it overrides nPeriods and stepsPerPeriod'),
        Field('stepsPerPeriod', 'Steps per period', 'int', '',
              'Integration steps per undulator period (ignored with a lattice file)'),
        Field('nPeriods', 'Periods', 'int', '',
              'Undulator periods (ignored with a lattice file)'),
        Field('sZ0', 'z₀', 'float', '', 'Starting zbar'),
        Field('iWriteNthSteps', 'Field dump every', 'int', 'steps',
              'Cadence for full field/electron dumps — these are the large files'),
        Field('iWriteIntNthSteps', 'Integrated dump every', 'int', 'steps',
              'Cadence for integrated data. Beware: a whole multiple of '
              'stepsPerPeriod samples every point on a period boundary, where '
              'the power dips sharply, and aliases with that dip'),
        Field('sRedistLen', 'Redistribute length', 'float', '',
              'Parallel field buffer length — how far ahead the buffer must reach'),
        Field('iRedistStp', 'Redistribute every', 'int', 'steps', ''),
        Field('iRandSeed', 'Random seed', 'int', '',
              'Fixes the shot-noise realisation. Leave at -1 to seed from the '
              'clock, which makes each run a different realisation'),
    ]),

    dict(name='Physics switches', file='main', block='mdata', fields=[
        Field('qOneD', '1D model', 'bool', '',
              'One field node and one macroparticle transversely'),
        Field('qFieldEvolve', 'Evolve field', 'bool', ''),
        Field('qElectronsEvolve', 'Evolve electrons', 'bool', ''),
        Field('qElectronFieldCoupling', 'Field → electrons', 'bool', ''),
        Field('qFocussing', 'Focusing', 'bool', ''),
        Field('qDiffraction', 'Diffraction', 'bool', ''),
        Field('qFilter', 'Filter', 'bool',
              '', 'Off means low frequencies are ignored during diffraction'),
        Field('qUndEnds', 'Undulator ends', 'bool', ''),
        Field('q_noise', 'Shot noise', 'bool',
              '', 'Shot noise in the initial electron distribution'),
        Field('qscaled', 'Scaled input', 'bool',
              '', 'Whether the values above are given scaled or in SI'),
    ]),

    dict(name='Beam', file='beam', block='blist', fields=[
        Field('Ipk', 'Peak current', 'float', 'A',
              'Sets the charge. Alternative to sQe'),
        Field('sQe', 'Bunch charge', 'float', 'C', 'Alternative to Ipk'),
        Field('sSigmaE', 'σ (x, y, z2, px, py, γ)', 'list', '',
              'Gaussian std dev in each dimension'),
        Field('sLenE', 'Length (x, y, z2, px, py, γ)', 'list', '',
              'Total modelled length in each dimension'),
        Field('iNumMPs', 'Macroparticles per dim', 'list', '',
              'Sampling in x, y, z2, px, py, γ. With the sequence loader only '
              'the z2 entry sets the slice count'),
        Field('nseqparts', 'Sequence particles', 'int', '',
              'Macroparticles per z2 slice, filling the other five dimensions'),
        Field('TrLdMeth', 'Loading method', 'choice:0|1|2', '',
              '0 equispaced, 1 random sequences, 2 Halton sequences'),
        Field('emitx', 'ε_x', 'float', 'm rad', 'Scaled transverse emittance'),
        Field('emity', 'ε_y', 'float', 'm rad', ''),
        Field('alphax', 'α_x', 'float', '', 'Twiss alpha'),
        Field('alphay', 'α_y', 'float', ''),
        Field('gammaf', 'γ / γ_r', 'float', '',
              'Ratio of average beam energy to the reference energy'),
        Field('chirp', 'Chirp', 'float', 'dγ/dz2', 'Energy chirp along z2'),
        Field('bcenter', 'Centre in z2', 'float', ''),
        Field('qMatched_A', 'Match to channel', 'bool', '',
              'Automatically match the beam to the focusing channel'),
    ]),

    dict(name='Seed', file='seed', block='slist', fields=[
        Field('sA0_X', 'A₀ x', 'float', '', 'Seed field amplitude, x'),
        Field('sA0_Y', 'A₀ y', 'float', '', 'Seed field amplitude, y'),
        Field('sSigmaF', 'σ (x, y, z2)', 'list', '', 'Seed field width'),
        Field('freqf', 'Frequency', 'float', '× ω_r', 'Seed frequency'),
        Field('meanZ2', 'Centre in z2', 'float', ''),
        Field('qFlatTop', 'Flat top', 'bool', '',
              'Flat-top rather than Gaussian longitudinal profile'),
    ]),
]


def all_fields():
    """(group, field) for every curated field, in display order."""
    for group in GROUPS:
        for field in group['fields']:
            yield group, field


def choices(kind):
    return kind.split(':', 1)[1].split('|') if kind.startswith('choice:') else None


# Labels for the choice fields whose raw values are opaque on their own.
CHOICE_LABELS = {
    'meshType': {'0': '0 — temporal', '1': '1 — periodic'},
    'TrLdMeth': {'0': '0 — equispaced', '1': '1 — random', '2': '2 — Halton'},
    'zundType': {'': '(1D — no off-axis variation)'},
}
