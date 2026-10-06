"""Radius sets stored in region definitions and checked by later steps (software tests).

Synthetic MD pair from ``test_cavity_obstacles``. ``region-prepare`` must write
schema version 2 with the radius set; ``region-trajectory``,
``hydration --cavity-results`` and ``region-compare`` must use that set or
refuse a different one unless explicitly overridden; version-1 definitions
must still load and mean the default set.
"""
import json

import numpy as np
import pytest

from crevice.cli import main
from crevice.radii import DEFAULT_RADII, radii_preset, reconcile_recorded_radii, same_radius_set
from crevice.region_definition import (declared_radii, load_definition, prepare_region, resolve_region_radii,
                                       trajectory_arguments, validate_definition, write_json)
from test_region_definition import definition


def prepared(tmp_path, **kwargs):
    path, data, _ = definition(tmp_path)
    files = prepare_region(path, tmp_path / 'review', **kwargs)
    return path, files


def test_prepare_writes_schema_2_with_the_default_set(tmp_path):
    _, files = prepared(tmp_path)
    data, _ = load_definition(files['definition_json'])
    assert data['schema_version'] == 2
    assert data['radii']['name'] == 'default' and same_radius_set(data['radii'], DEFAULT_RADII)
    manifest = json.loads((tmp_path / 'review' / 'manifest.json').read_text())
    assert manifest['radii']['name'] == 'default'
    assert manifest['radii_check']['status'] == 'no_declared_set'
    review = json.loads((tmp_path / 'review' / 'synthetic_review.json').read_text())
    assert review['radii']['table_sha256'] == DEFAULT_RADII.table_sha256


def test_prepare_with_an_explicit_set_stores_it_and_trajectory_uses_it(tmp_path):
    _, files = prepared(tmp_path, radii='hole')
    data, _ = load_definition(files['definition_json'])
    assert data['radii']['name'] == 'hole'
    argv, metadata = trajectory_arguments(files['definition_json'], tmp_path / 'run')
    assert metadata['radius_set'].name == 'hole'
    assert metadata['radii_check']['status'] == 'definition_set_used'
    out = tmp_path / 'run'
    assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out), '--skip-hydration',
                 '--stop', '2', '--dpi', '60']) == 0
    report = json.loads((out / 'synthetic_cavity_statistics.json').read_text())
    assert report['settings']['radii']['name'] == 'hole'
    assert report['settings']['region_definition']['radii']['name'] == 'hole'
    manifest = json.loads((out / 'synthetic_cavity_manifest.json').read_text())
    assert manifest['radii']['name'] == 'hole'


def test_region_trajectory_refuses_a_different_set_unless_overridden(tmp_path, capsys):
    _, files = prepared(tmp_path)
    out = tmp_path / 'run'
    assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out), '--skip-hydration',
                 '--stop', '2', '--radii', 'hole']) == 2
    assert 'Radius set mismatch' in capsys.readouterr().err
    assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out), '--skip-hydration',
                 '--stop', '2', '--dpi', '60', '--radii', 'hole', '--allow-radii-mismatch']) == 0
    provenance = json.loads((out / 'synthetic_region_provenance.json').read_text())
    assert provenance['radii']['name'] == 'hole'
    assert provenance['radii_check']['status'] == 'mismatch_overridden'
    assert provenance['radii_check']['declared']['name'] == 'default'
    # The same set given explicitly is accepted.
    out2 = tmp_path / 'run2'
    assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out2), '--skip-hydration',
                 '--stop', '2', '--dpi', '60', '--radii', 'default']) == 0


def test_version_1_definitions_still_load_and_mean_the_default_set(tmp_path):
    _, files = prepared(tmp_path)
    path = files['definition_json']
    data = json.loads(open(path).read())
    data['schema_version'] = 1
    data.pop('radii')
    write_json(path, data)
    loaded, _ = load_definition(path)
    record, origin = declared_radii(loaded)
    assert origin == 'legacy_v1_default' and same_radius_set(record, DEFAULT_RADII)
    used, check = resolve_region_radii(loaded)
    assert used.table_sha256 == DEFAULT_RADII.table_sha256 and check['status'] == 'definition_set_used'
    with pytest.raises(ValueError, match='Radius set mismatch'):
        resolve_region_radii(loaded, 'bondi')
    # Unprepared version-1 requests have not chosen a set: any requested set is used.
    request, _, _ = definition(tmp_path / 'fresh')
    raw, _ = load_definition(request)
    assert declared_radii(raw) == (None, 'not_declared')
    assert resolve_region_radii(raw, 'hole')[0].name == 'hole'


@pytest.mark.parametrize('change, match', [
    (lambda d: d.pop('radii'), 'missing fields'),
    (lambda d: d.update(radii={'name': 'x'}), 'table_sha256'),
    (lambda d: d.update(schema_version=3), 'schema_version'),
])
def test_schema_2_requires_a_radius_record(tmp_path, change, match):
    _, files = prepared(tmp_path)
    data = json.loads(open(files['definition_json']).read())
    change(data)
    with pytest.raises(ValueError, match=match):
        validate_definition(data)
    legacy = json.loads(open(files['definition_json']).read())
    legacy['schema_version'] = 1
    with pytest.raises(ValueError, match='unknown fields'):
        validate_definition(legacy)


def test_hydration_cavity_results_uses_or_checks_the_recorded_set(tmp_path):
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    _, files = prepared(tmp_path, radii='hole')
    out = tmp_path / 'run'
    assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out), '--skip-hydration',
                 '--stop', '2', '--dpi', '60']) == 0
    stats = out / 'synthetic_cavity_statistics.json'
    reader = json.loads((out / 'synthetic_reader.json').read_text())
    with pytest.raises(ValueError, match='Radius set mismatch'):
        analyze_hydration_trajectory(reader['topology'], reader['trajectory'], context_json=stats, radii='bondi')
    record = json.loads(stats.read_text())['settings']['radii']
    used, check = reconcile_recorded_radii(record, None)
    assert used.name == 'hole' and check['status'] == 'recorded_set_used'
    used, check = reconcile_recorded_radii(record, 'bondi', allow_mismatch=True)
    assert used.name == 'bondi' and check['status'] == 'mismatch_overridden'
    used, check = reconcile_recorded_radii(None, None)
    assert used is DEFAULT_RADII and check['status'] == 'unverified_no_recorded_set'


def test_region_compare_refuses_runs_with_different_radius_sets(tmp_path):
    from crevice.region_comparison import compare_regions
    path, data, _ = definition(tmp_path)
    stats = []
    for name, radii in [('left', 'default'), ('right', 'hole')]:
        data['region_id'] = name
        write_json(path, data)
        files = prepare_region(path, tmp_path / (name + '-review'), radii=radii)
        out = tmp_path / name
        assert main(['region-trajectory', files['definition_json'], '--out-dir', str(out), '--skip-hydration',
                     '--stop', '2', '--dpi', '60']) == 0
        stats.append(out / (name + '_cavity_statistics.json'))
    with pytest.raises(ValueError, match='Radius set mismatch'):
        compare_regions(stats[0], stats[1], tmp_path / 'cmp')
    result = compare_regions(stats[0], stats[1], tmp_path / 'cmp2', allow_radii_mismatch=True)
    assert result['radii']['status'] == 'mismatch_overridden'
    assert result['radii']['right']['name'] == 'hole'
    left = np.load(tmp_path / 'left' / 'left_cavity_frames.npz')['volume_A3']
    right = np.load(tmp_path / 'right' / 'right_cavity_frames.npz')['volume_A3']
    assert np.all(right < left)  # larger HOLE carbon radii give smaller regional volumes
