"""Figures and viewer scenes are clean by default; one opt-in restores their text.

Software tests on synthetic structures: default figures carry no titles,
suptitles or in-plot notes, default scene scripts contain no text commands,
``annotate=True``/``--annotate`` restores them, the text is recorded in PNG
metadata and scene JSON either way, and every data output is unchanged.
"""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
import json
import re
from pathlib import Path

import numpy as np
import pytest

import crevice.figures as figures
from crevice.cli import build_parser, main
from crevice.presentation import (ANNOTATION_METADATA_KEY, annotations_enabled, figure_annotations,
                                  figure_description, write_viewer_structure)
from test_channel_coordinate import rectangular_channel

# Category labels drawn as text artists: they are the chart's tick labels.
CATEGORY_LABEL_FIGURES = ('_network_chord.png', '_connectivity.png', '_analysis_overview.png')

# Commands that draw text in native viewer scripts.
SCENE_TEXT = [re.compile(p, re.M) for p in (
    r"pseudoatom\([^\n]*label=",       # PyMOL label pseudoatoms
    r"show\('labels'",                  # PyMOL label representation
    r"\bcmd\.label\(",
    r"\b2dlabels\b",                    # ChimeraX 2D labels
    r"run\(session,'label ",            # ChimeraX 3D atom labels
    r"run\(session,'key ",              # ChimeraX colour key
    r"^\s*key ",
    r"graphics \S+ text ",              # VMD text graphics
)]
SCENE_SUFFIXES = ('.pml', '.vmd', '.tcl', '.cxc', '_chimerax.py')


def _figure_text(fig):
    titles = [t.get_text() for ax in fig.axes for t in (ax.title, ax._left_title, ax._right_title) if t.get_text()]
    return {'titles': titles,
            'suptitle': fig._suptitle.get_text() if fig._suptitle is not None else '',
            'texts': [t.get_text() for ax in fig.axes for t in ax.texts if t.get_text()],
            'figure_texts': [t.get_text() for t in fig.texts if t.get_text()]}


@pytest.fixture
def saved_figures(monkeypatch):
    """Record the text artists of every figure CREVICE saves."""
    seen = {}
    original = figures._save_figure

    def spy(fig, path, *, dpi, text=None):
        seen[Path(path).name] = _figure_text(fig)
        return original(fig, path, dpi=dpi, text=text)
    monkeypatch.setattr(figures, '_save_figure', spy)
    return seen


def _notices(png):
    return {line.split(': ', 1)[1] for line in figure_description(png)['Description'].splitlines() if line.startswith('Notice: ')}


def _assert_clean(name, text, png):
    assert not text['titles'] and not text['suptitle'] and not text['figure_texts'], (name, text)
    if not name.endswith(CATEGORY_LABEL_FIGURES):
        allowed = {' '.join(n.split()) for n in _notices(png)}
        assert {' '.join(t.split()) for t in text['texts']} <= allowed, (name, text['texts'])
    assert figure_description(png)[ANNOTATION_METADATA_KEY].startswith('omitted')


def _normalised(path, root):
    return path.read_bytes().replace(str(root).encode(), b'ROOT')


def _publish(tmp_path, name, **options):
    frame = rectangular_channel()
    source = write_viewer_structure(frame, tmp_path/'channel.pdb')
    root = tmp_path/name
    figures.write_static_publication_bundle(frame, structure_path=source, output_dir=root, prefix='c', samples=41, enclosure_radius=0.8,
                                            search_radius=8, cast_spacing=.5, dpi=60,
                                            residue_groups={a.residue_key.label: 'H1' if a.x < 0 else 'H2' for a in frame.atoms},
                                            **options)
    return root


def test_publish_figures_are_clean_by_default_and_annotate_restores_text(tmp_path, saved_figures):
    clean = _publish(tmp_path, 'clean')
    default = dict(saved_figures); saved_figures.clear()
    drawn = _publish(tmp_path, 'annotated', annotate=True)
    annotated = dict(saved_figures)
    assert set(default) == set(annotated) and len(default) >= 7
    for name, text in default.items():
        _assert_clean(name, text, clean/name)
        # The omitted text is recorded in the file itself.
        description = figure_description(clean/name)
        assert description['Title'] and description['Description'] and description['Software'] == 'CREVICE'
        if not _notices(clean/name):   # an empty-data panel never had a title
            assert annotated[name]['titles'] or annotated[name]['suptitle'], name
        assert figure_description(drawn/name)[ANNOTATION_METADATA_KEY] == 'drawn'
    # Value call-outs, mouth labels and residue landmark names return with the opt-in.
    assert any(t.endswith(' Å') for t in annotated['c_profile_radius.png']['texts'])
    landmarks = [line for line in figure_description(clean/'c_profile_radius_annotated.png')['Description'].splitlines()
                 if line.startswith('Note: Residue landmarks')]
    assert landmarks and 'A:ALA' in landmarks[0]
    assert any(t.startswith('A:ALA') for t in annotated['c_profile_radius_annotated.png']['texts'])
    assert not any(t.startswith('A:ALA') for t in default['c_profile_radius_annotated.png']['texts'])
    # Chord node names are category labels: drawn either way.
    assert default['c_network_chord.png']['texts'] == annotated['c_network_chord.png']['texts']


def test_publish_data_outputs_identical_with_and_without_annotations(tmp_path):
    clean, drawn = _publish(tmp_path, 'clean'), _publish(tmp_path, 'annotated', annotate=True)
    names = sorted(p.relative_to(clean) for p in clean.rglob('*') if p.is_file())
    assert names == sorted(p.relative_to(drawn) for p in drawn.rglob('*') if p.is_file())
    compared = [n for n in names if n.suffix != '.png']
    assert len(compared) > 20
    for name in compared:
        assert _normalised(clean/name, clean) == _normalised(drawn/name, drawn), name


def _toy_frame():
    from crevice import Atom, StructureFrame
    def atom(i, name, xyz, *, resid=1, resname='ALA', element='N', hetero=False):
        return Atom(i, name, resname, 'A', resid, *xyz, element=element, hetero=hetero)
    return StructureFrame((atom(1, 'N', (0, 0, 0)), atom(2, 'CA', (1.45, 0, 0), element='C'),
                           atom(3, 'C', (2, 1.4, 0), element='C'), atom(4, 'O', (1.8, 2.5, 0), element='O'),
                           atom(5, 'O', (0, 0, 3), resid=2, resname='HOH', element='O', hetero=True)))


def _scene_files(root):
    return [p for p in root.iterdir() if p.name.endswith(SCENE_SUFFIXES)]


def _text_commands(path):
    source = path.read_text()
    return [p.pattern for p in SCENE_TEXT if p.search(source)]


def test_hydration_scenes_have_no_text_by_default_and_annotate_restores_it(tmp_path, saved_figures):
    from crevice.hydration_workflow import static_hydration_bundle
    clean, drawn = tmp_path/'clean', tmp_path/'annotated'
    static_hydration_bundle(_toy_frame(), clean, prefix='toy', sasa_samples=32, dpi=60)
    default = dict(saved_figures); saved_figures.clear()
    static_hydration_bundle(_toy_frame(), drawn, prefix='toy', sasa_samples=32, dpi=60, annotate=True)
    for name, text in default.items():
        _assert_clean(name, text, clean/name)
    scripts = _scene_files(clean)
    assert len(scripts) >= 8
    for path in scripts:
        assert not _text_commands(path), (path.name, _text_commands(path))
    annotated = {p.name: _text_commands(p) for p in _scene_files(drawn)}
    assert annotated['toy_analysis.pml'] and annotated['toy_analysis_chimerax.py'] and annotated['toy_analysis.vmd']
    scene = json.loads((clean/'toy_analysis_scene.json').read_text())
    assert scene['annotations']['drawn'] is False
    assert scene['annotations']['texts']['colour_scale_caption'] == 'Observed snapshot contact'
    assert json.loads((drawn/'toy_analysis_scene.json').read_text())['annotations']['drawn'] is True
    # No HTML report is written; figure descriptions stay in the PNG metadata.
    assert not list(clean.glob('*.html')) and not list(drawn.glob('*.html'))
    assert 'water-contact occupancy' in figure_description(clean/'toy_analysis_overview.png')['Description']
    # Data outputs are unchanged; only images, scene text and their JSON records differ.
    presentation = {'toy_analysis_scene.json'}
    for path in sorted(clean.iterdir()):
        if path.suffix == '.png' or path.name.endswith(SCENE_SUFFIXES) or path.name in presentation:
            continue
        assert _normalised(path, clean) == _normalised(drawn/path.name, drawn), path.name
    other = json.loads((drawn/'toy_analysis_scene.json').read_text())
    scene.pop('annotations'); other.pop('annotations')
    assert json.dumps(scene).replace(str(clean), '') == json.dumps(other).replace(str(drawn), '')


def test_region_review_scene_labels_are_opt_in(tmp_path):
    from crevice.region_definition import prepare_region
    from test_region_definition import definition
    path, _, _ = definition(tmp_path)
    clean = prepare_region(path, tmp_path/'clean')
    drawn = prepare_region(path, tmp_path/'annotated', annotate=True)
    for key in ['review_pml', 'review_vmd', 'review_tcl', 'review_chimerax_python']:
        assert not _text_commands(Path(clean[key])), key
        assert _text_commands(Path(drawn[key])), key
    scene = json.loads(Path(clean['review_scene_json']).read_text())
    assert scene['annotations']['drawn'] is False and scene['annotations']['texts']['residue_labels']
    assert 'review_html' not in clean and not list((tmp_path/'clean').glob('*.html'))
    assert Path(clean['review_landmarks_csv']).read_text().splitlines()[0].startswith('residue,anchor_atoms,nearest_reference_sample_A')


def test_member_water_group_has_no_text(tmp_path):
    from crevice.analysis_scene_templates import write_scene_scripts, CXC_PY, PML
    for template in (PML, CXC_PY):
        body = template.replace('__LEGEND__', '').replace('__REGION_STATUS__', '').replace('__LABEL_VISIBILITY__', '')
        assert 'crevice_member_waters' in body
        assert not [p.pattern for p in SCENE_TEXT if p.search(body)]


def test_explicit_argument_overrides_context_and_caller_title_is_always_drawn(tmp_path, saved_figures):
    from crevice import pore_profile
    profile = pore_profile(rectangular_channel(), enclosure_radius=0.8, samples=21, search_radius=8)
    with figure_annotations():
        assert annotations_enabled()
        figures.plot_profile_radius(profile, tmp_path/'forced_clean.png', dpi=40, annotate=False)
        figures.plot_profile_radius(profile, tmp_path/'inherited.png', dpi=40)
    assert not annotations_enabled()
    assert figure_description(tmp_path/'forced_clean.png')[ANNOTATION_METADATA_KEY].startswith('omitted')
    assert figure_description(tmp_path/'inherited.png')[ANNOTATION_METADATA_KEY] == 'drawn'
    figures.plot_profile_radius(profile, tmp_path/'titled.png', dpi=40, title='My channel')
    assert saved_figures['titled.png']['titles'] == ['My channel']
    assert saved_figures['forced_clean.png']['titles'] == [] and saved_figures['inherited.png']['titles']


def test_cli_annotate_is_one_flag_on_every_figure_command(tmp_path):
    parser = build_parser()
    commands = parser._subparsers._group_actions[0].choices
    figure_commands = {'region-prepare', 'region-compare', 'region-trajectory', 'hydration', 'cavity-trajectory',
                       'cast', 'profile', 'residues', 'cavities', 'network', 'residue-evidence',
                       'trajectory', 'publish', 'analyze', 'static-suite'}
    for name, sub in commands.items():
        flags = {s for action in sub._actions for s in action.option_strings}
        assert ('--annotate' in flags) == (name in figure_commands), name
    assert '--annotate' in parser.format_help()
    source = write_viewer_structure(rectangular_channel(), tmp_path/'channel.pdb')
    for name, extra in [('clean', []), ('annotated', ['--annotate'])]:
        assert main(['profile', str(source), '--offline', '--enclosure-radius', '0.8', '--samples', '21', '--search-radius', '8', '-o',
                     str(tmp_path/(name+'.json')), '--png', str(tmp_path/(name+'.png')),
                     '--annotated-png', str(tmp_path/(name+'_landmarks.png')), '--dpi', '40', *extra]) == 0
    assert not annotations_enabled()
    assert figure_description(tmp_path/'clean.png')[ANNOTATION_METADATA_KEY].startswith('omitted')
    assert figure_description(tmp_path/'annotated_landmarks.png')[ANNOTATION_METADATA_KEY] == 'drawn'
    assert (tmp_path/'clean.json').read_bytes() == (tmp_path/'annotated.json').read_bytes()
