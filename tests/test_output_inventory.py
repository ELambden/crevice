"""Default output inventories: CSV tables, no HTML, no cast preview or tunnel summary.

Software tests on synthetic channels only; no biological result is implied.
They check which files each command writes by default, that every tabular
dataset has a CSV whose column names carry units, that the JSON files that
remain are the documented structured records, and that the residue-contact and
network figures explain their colours.
"""
import csv
import json
from pathlib import Path

import pytest

import crevice.figures as figures
from crevice.cli import main
from crevice.presentation import write_viewer_structure
from test_channel_coordinate import rectangular_channel

# Never written by default any more.
FORBIDDEN_SUFFIXES = ('.html', '_cast_preview.png', '_tunnel_summary.png', '_analysis_report_data.json',
                      '_residue_contacts.json', '_residues.json', '_features.json', '_hydration_mobility.json',
                      '_profile_sensitivity.json', '_tunnels.png')
# JSON kept as structured records, provenance or viewer scene definitions (see docs/getting-started/outputs.md).
KEPT_JSON_SUFFIXES = ('_manifest.json', '_input_report.json', '_profile.json', '_profile_status.json',
                      '_void_cast.json', '_selection.json', '_region_annotations.json', '_cavities.json',
                      '_tunnels.json', '_network.json', '_connectivity.json', '_scene.json',
                      '_volume_metadata.json', '.identities.json', '_hydration.json', '_interactions.json',
                      '_water_density.json', '_residue_evidence.json', '_hydration_alignment.json')
PROFILE_COLUMNS = ['index', 'x_A', 'y_A', 'z_A', 'radius_A', 'raw_clearance_A', 'axis_position_A',
                   'nearest_atom_serial', 'nearest_residue']


def _header(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return next(csv.reader(handle))


def _rows(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def _all_files(root):
    return sorted(p for p in Path(root).rglob('*') if p.is_file())


def _assert_clean_inventory(root):
    files = _all_files(root)
    assert files
    for path in files:
        assert not path.name.endswith(FORBIDDEN_SUFFIXES), path.name
        if path.suffix == '.json':
            assert path.name.endswith(KEPT_JSON_SUFFIXES), path.name


@pytest.fixture(scope='module')
def channel_pdb(tmp_path_factory):
    root = tmp_path_factory.mktemp('channel')
    return write_viewer_structure(rectangular_channel(), root / 'channel.pdb')


PROFILE_ARGS = ['--offline', '--samples', '21', '--search-radius', '8']


@pytest.fixture(scope='module')
def published(tmp_path_factory, channel_pdb):
    root = tmp_path_factory.mktemp('publish')
    assert main(['publish', str(channel_pdb), *PROFILE_ARGS, '--cast-spacing', '.5', '--dpi', '40',
                 '--hydration-sasa-points', '32', '--out-dir', str(root), '--prefix', 'c']) == 0
    return root


def test_publish_default_inventory_is_csv_first_without_html(published):
    _assert_clean_inventory(published)
    manifest = json.loads((published / 'c_manifest.json').read_text())
    files = manifest['files']
    for key in ['profile_csv', 'residue_contacts_csv', 'void_cast_csv', 'cavities_csv', 'tunnels_csv',
                'tunnel_points_csv', 'network_nodes_csv', 'network_edges_csv', 'hydration_residues_csv']:
        assert Path(files[key]).is_file(), key
    assert 'tunnel_summary_png' not in files and 'residue_contacts_json' not in files
    assert _header(files['profile_csv']) == PROFILE_COLUMNS
    assert _header(files['residue_contacts_csv'])[:3] == ['residue', 'min_distance_A', 'mean_distance_A']
    assert {'volume_A3', 'center_x_A', 'max_radius_A'} <= set(_header(files['void_cast_csv']))
    assert {'cavity_id', 'volume_A3', 'max_radius_A'} <= set(_header(files['cavities_csv']))
    assert {'tunnel_id', 'length_A', 'bottleneck_radius_A'} <= set(_header(files['tunnels_csv']))
    assert _header(files['tunnel_points_csv']) == ['tunnel_id', *PROFILE_COLUMNS]
    assert {'node_id', 'kind', 'degree', 'residue_degree', 'x_A'} <= set(_header(files['network_nodes_csv']))
    assert {'source', 'target', 'interaction', 'distance_A', 'cutoff_A'} <= set(_header(files['network_edges_csv']))
    # The CSVs hold the same numbers as the retained JSON records.
    profile = json.loads(Path(files['profile_json']).read_text())
    assert min(float(r['radius_A']) for r in _rows(files['profile_csv'])) == pytest.approx(profile['min_radius'], abs=0)
    volumes = [float(r['volume_A3']) for r in _rows(files['void_cast_csv'])]
    assert sum(volumes) == pytest.approx(manifest['void_cast']['total_volume'], abs=1e-9)
    tunnels = json.loads(Path(files['tunnels_json']).read_text())['tunnels']
    assert len(_rows(files['tunnels_csv'])) == len(tunnels)
    assert len(_rows(files['tunnel_points_csv'])) == sum(len(t['points']) for t in tunnels)
    network = json.loads(Path(files['network_json']).read_text())
    assert len(_rows(files['network_edges_csv'])) == network['network']['edge_count']
    degrees = {r['node_id']: int(r['degree']) for r in _rows(files['network_nodes_csv'])}
    assert degrees == network['metrics']['degree']


def test_hydration_views_write_frame_tables_and_no_html(published):
    hydration = json.loads((published / 'c_hydration_manifest.json').read_text())['files']
    assert _header(hydration['hydration_frames_csv'])[:3] == ['frame_index', 'time_ps', 'residue']
    assert {'water_count', 'sasa_A2'} <= set(_header(hydration['hydration_frames_csv']))
    assert _header(hydration['hydration_region_frames_csv'])[:2] == ['frame_index', 'time_ps']
    assert _header(hydration['interaction_edge_frames_csv']) == ['frame_row', 'time_ps', 'kind', 'source', 'target',
                                                                 'observation_count']
    assert {'phi_degrees', 'chi1_degrees', 'secondary_structure'} <= set(_header(hydration['residue_conformation_frames_csv']))
    assert _header(hydration['hydration_shared_water_csv'])[:2] == ['source', 'target']
    assert 'analysis_html' not in hydration and not list(published.glob('*.html'))


def test_analyze_single_commands_and_trajectory_write_csv_companions(tmp_path, channel_pdb):
    out = tmp_path / 'analyze'
    assert main(['analyze', str(channel_pdb), *PROFILE_ARGS, '--analysis', 'profile', '--analysis', 'residues',
                 '--analysis', 'cavities', '--analysis', 'tunnels', '--analysis', 'network', '--analysis', 'features',
                 '--png', '--dpi', '40', '--skip-hydration', '--out-dir', str(out), '--prefix', 'c']) == 0
    _assert_clean_inventory(out)
    files = json.loads((out / 'c_manifest.json').read_text())['files']
    assert {'profile.json', 'profile.csv', 'residues.csv', 'cavities.json', 'cavities.csv', 'tunnels.json',
            'tunnels.csv', 'tunnel_points.csv', 'network.json', 'network_nodes.csv', 'network_edges.csv',
            'features.csv', 'profile.png', 'residues.png', 'cavities.png', 'network.png'} == set(files)
    assert _header(files['features.csv']) == ['feature', 'value', 'unit']
    units = {r['feature']: r['unit'] for r in _rows(files['features.csv'])}
    assert units['min_radius'] == 'Å' and units['volume_estimate'] == 'Å³' and units['contact_count'] == ''

    single = tmp_path / 'single'
    assert main(['cavities', str(channel_pdb), '--offline', '-o', str(single / 'cav.json')]) == 0
    assert main(['tunnels', str(channel_pdb), '--offline', '-o', str(single / 'tun.json')]) == 0
    assert main(['network', str(channel_pdb), *PROFILE_ARGS, '-o', str(single / 'net.json'),
                 '--connectivity-json', str(single / 'con.json')]) == 0
    assert main(['features', str(channel_pdb), *PROFILE_ARGS, '-o', str(single / 'feat.json')]) == 0
    assert sorted(p.name for p in single.iterdir()) == sorted([
        'cav.json', 'cav.csv', 'tun.json', 'tun.csv', 'tun_points.csv', 'net.json', 'net_nodes.csv',
        'net_edges.csv', 'con.json', 'con.csv', 'con_groups.csv', 'feat.json', 'feat.csv'])
    assert _header(single / 'con.csv')[:3] == ['residue', 'group', 'is_lining']
    with pytest.raises(SystemExit):   # the tunnel summary figure was removed
        main(['tunnels', str(channel_pdb), '--offline', '-o', str(single / 'x.json'), '--png', str(single / 'x.png')])

    trajectory = tmp_path / 'trajectory'
    assert main(['trajectory', str(channel_pdb), str(channel_pdb), '-o', str(trajectory / 't.json'),
                 '--samples', '21', '--search-radius', '8', '--no-align', '--confidence', '0',
                 '--distribution-csv', str(trajectory / 'd.csv'), '--skip-hydration']) == 0
    frames = _rows(trajectory / 't_frames.csv')
    assert [r['frame_index'] for r in frames] == ['0', '1'] and frames[0]['min_radius_A']
    assert _header(trajectory / 't_profiles.csv')[:3] == ['frame_index', 'time_ps', 'index']
    assert len(_rows(trajectory / 't_profiles.csv')) == 2 * int(frames[0]['profile_point_count'])
    assert _header(trajectory / 'd.csv')[0] == 'coordinate_A' and 'mean_radius_A' in _header(trajectory / 'd.csv')


def test_cast_writes_region_table_and_no_preview(tmp_path):
    source = write_viewer_structure(rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True),
                                    tmp_path / 'pocket.pdb')
    out = tmp_path / 'cast'
    assert main(['cast', str(source), '--offline', '--out-dir', str(out), '--prefix', 'p', '--min-component-volume', '30',
                 '--skip-hydration']) == 0
    _assert_clean_inventory(out)
    manifest = json.loads((out / 'p_manifest.json').read_text())
    assert 'cast_preview_png' not in manifest['files'] and not list(out.glob('*.png'))
    rows = _rows(manifest['files']['void_cast_csv'])
    assert rows and sum(float(r['volume_A3']) for r in rows) == pytest.approx(manifest['void_cast']['total_volume'], abs=1e-9)


# --- figures ------------------------------------------------------------------

def _snapshot(fig):
    """What a reader sees: legends, axis labels, ticks and bar geometry (the figure is cleared on save)."""
    from matplotlib.colors import to_hex
    def legend(item):
        return {'labels': [t.get_text() for t in item.get_texts()],
                'colours': [to_hex(h.get_facecolor() if hasattr(h, 'get_facecolor') else h.get_color())
                            for h in item.legend_handles]}
    axes = []
    for ax in fig.axes:
        axes.append({'xlabel': ax.get_xlabel(), 'yticks': [t.get_text() for t in ax.get_yticklabels()],
                     'bars': [(p.get_x(), p.get_y() + p.get_height() / 2, p.get_width(), to_hex(p.get_facecolor()))
                              for p in ax.patches if hasattr(p, 'get_x') and hasattr(p, 'get_width')],
                     'legend': legend(ax.get_legend()) if ax.get_legend() else None})
    return {'legends': [legend(item) for item in fig.legends], 'axes': axes}


@pytest.fixture
def saved(monkeypatch):
    seen = {}
    original = figures._save_figure

    def spy(fig, path, *, dpi, text=None):
        seen[Path(path).name] = _snapshot(fig)
        return original(fig, path, dpi=dpi, text=text)
    monkeypatch.setattr(figures, '_save_figure', spy)
    return seen


def test_residue_contacts_legend_names_exactly_the_role_colours(tmp_path, saved):
    from crevice.models import ResidueContact
    roles = ['bottleneck', 'bottleneck-nearby', 'lining', 'nearby']
    contacts = [ResidueContact(f'A:ALA{i}', 0.5 * i, 1.0, 5, 3, role, ('hydrophobic',), 5.0 - i)
                for i, role in enumerate(roles)]
    figures.plot_residue_contacts(contacts, tmp_path / 'contacts.png', dpi=40)
    fig = saved['contacts.png']
    assert len(fig['legends']) == 1
    legend = fig['legends'][0]
    assert legend['labels'] == [label for _role, label in figures.ROLE_LEGEND]
    assert legend['colours'] == [figures._role_color(role).lower() for role, _label in figures.ROLE_LEGEND]
    # Every bar's colour is one the legend explains.
    assert {bar[3] for bar in fig['axes'][0]['bars']} <= set(legend['colours'])
    assert fig['axes'][0]['xlabel'].startswith('Influence score (')
    # Only roles that are drawn appear in the legend.
    figures.plot_residue_contacts(contacts[2:3], tmp_path / 'one.png', dpi=40)
    assert saved['one.png']['legends'][0]['labels'] == ['Lining (gap ≤ 1 Å elsewhere)']


def test_network_summary_axis_is_explained_and_counts_residue_contacts(tmp_path, saved):
    from crevice import annotate_residues, build_cavity_network, pore_profile
    from crevice.io import write_network_csv
    frame = rectangular_channel()
    profile = pore_profile(frame, samples=21, search_radius=8)
    network = build_cavity_network(frame, profile, contacts=annotate_residues(frame, profile.points))
    figures.plot_network_summary(network, tmp_path / 'summary.png', dpi=40)
    fig = saved['summary.png']
    ax = fig['axes'][0]
    assert ax['xlabel'] == 'Residue contacts (degree: residues with atom-surface gap ≤ 4.5 Å)'
    legend = fig['legends'][0]
    assert legend['labels'] and set(legend['labels']) <= {label for _k, label in figures.INTERACTION_LEGEND}
    # Stacked bar lengths per residue equal its residue-residue degree in the node table.
    write_network_csv(network, tmp_path / 'nodes.csv', tmp_path / 'edges.csv')
    degree = {r['node_id']: int(r['residue_degree']) for r in _rows(tmp_path / 'nodes.csv')}
    totals = {}
    for x, y, width, _colour in ax['bars']:
        label = ax['yticks'][round(y)]
        totals[label] = max(totals.get(label, 0), x + width)
    assert totals and all(totals[k] == degree[k] for k in totals)
    assert not any(k.startswith('region:') for k in ax['yticks'])
    # Bar colours are explained by the legend.
    assert {bar[3] for bar in ax['bars']} <= set(legend['colours'])


def test_chord_colours_follow_the_interaction_names_networks_write(tmp_path, saved):
    written = ['salt_bridge_candidate', 'oppositely_charged_residue_contact', 'hydrophobic_residue_contact',
               'polar_residue_contact', 'distance_contact', 'region_lining', 'region_bottleneck']
    colours = {name: figures._interaction_color(name) for name in written}
    assert colours['region_lining'] == colours['region_bottleneck'] == figures.CREVICE_RED
    assert len({colours[n] for n in written[:5]}) == 5   # each residue class has its own colour
    from crevice import annotate_residues, build_cavity_network, pore_profile
    frame = rectangular_channel()
    profile = pore_profile(frame, samples=21, search_radius=8)
    network = build_cavity_network(frame, profile, contacts=annotate_residues(frame, profile.points))
    figures.plot_network_chord(network, tmp_path / 'chord.png', dpi=40)
    # The chord legend is a figure legend in its own band below the circle.
    assert saved['chord.png']['axes'][0]['legend'] is None
    labels = saved['chord.png']['legends'][0]['labels']
    drawn = {figures._interaction_key(e.interaction) for e in network.edges}
    assert set(labels) <= {label for key, label in figures.INTERACTION_LEGEND if key in drawn}
    assert 'Channel/cavity contact' in labels
