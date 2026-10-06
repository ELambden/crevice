"""Nested interval statistics of the trajectory summaries are also written as flat CSV."""
import csv
import json
from pathlib import Path

from crevice.io import SUMMARY_STATISTICS_FIELDS, summary_statistics_rows, write_summary_statistics_csv


def read_rows(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def test_flattening_keeps_every_interval_record_and_names_list_items(tmp_path):
    stat = lambda mean, lower=None: {  # noqa: E731
        'total_frames': 4, 'observed_frames': 4, 'missing_frames': 0, 'coverage_fraction': 1.0, 'mean': mean,
        'confidence_interval': {'status': 'complete' if lower is not None else 'unavailable_observed_constant',
                                'method': 'joint_studentized_batch_bootstrap', 'confidence': .95,
                                'pointwise_lower': [lower], 'pointwise_upper': [None if lower is None else mean + 1]}}
    report = {'settings': {'probe_radius': 1.4}, 'volume_A3': stat(100., 90.),
              'residues': [{'residue': 'A:ALA1', 'boundary_area_A2': stat(5., 4.), 'volume_correlation': None},
                           {'residue': 'A:GLY2', 'boundary_area_A2': stat(0.)}],
              'contacts': [{'source': 'A:ALA1', 'target': 'A:GLY2', 'occupancy_statistics': stat(.5, .25)}],
              'correlations': {'volume_vs_count': .75, 'interpretation': 'text is not a row'}}
    rows = summary_statistics_rows(report)
    assert [(r['quantity'], r['item']) for r in rows] == [
        ('volume_A3', ''), ('residues.boundary_area_A2', 'A:ALA1'), ('residues.volume_correlation', 'A:ALA1'),
        ('residues.boundary_area_A2', 'A:GLY2'), ('contacts.occupancy_statistics', 'A:ALA1|A:GLY2'),
        ('correlations.volume_vs_count', '')]
    assert rows[0]['mean_lower'] == 90. and rows[0]['mean_upper'] == 101. and rows[0]['interval_status'] == 'complete'
    assert rows[3]['mean_lower'] == '' and rows[3]['interval_status'] == 'unavailable_observed_constant'
    assert rows[2]['value'] == '' and rows[5]['value'] == .75 and rows[5]['mean'] == ''
    write_summary_statistics_csv(report, tmp_path / 'summary.csv')
    fields, written = read_rows(tmp_path / 'summary.csv')
    assert tuple(fields) == SUMMARY_STATISTICS_FIELDS and len(written) == 6
    assert float(written[4]['mean_upper']) == 1.5
    assert write_summary_statistics_csv({'nothing': 1}, tmp_path / 'empty.csv') == []
    fields, written = read_rows(tmp_path / 'empty.csv')
    assert tuple(fields) == SUMMARY_STATISTICS_FIELDS and written == []


def test_every_json_interval_statistic_appears_once_in_the_bundle_csvs(tmp_path):
    """Cavity-trajectory, hydration and region-comparison bundles all write the flat table."""
    from test_region_definition import definition
    from crevice.cli import main
    from crevice.region_comparison import compare_regions
    from crevice.region_definition import prepare_region, write_json

    path, data, _ = definition(tmp_path)
    runs = []
    for name in ['left', 'right']:
        data['region_id'] = name
        write_json(path, data)
        prepared = prepare_region(path, tmp_path / (name + '-review'))
        out = tmp_path / name
        assert main(['region-trajectory', prepared['definition_json'], '--out-dir', str(out), '--stop', '2', '--dpi', '60']) == 0
        runs.append(out)

    def count_statistics(node):
        if isinstance(node, dict):
            if 'mean' in node and 'confidence_interval' in node:
                return 1
            return sum(count_statistics(v) for v in node.values())
        if isinstance(node, list):
            return sum(count_statistics(v) for v in node)
        return 0

    for json_name, csv_name in [('left_cavity_statistics.json', 'left_cavity_statistics_summary.csv'),
                                ('left_hydration.json', 'left_hydration_summary.csv')]:
        report = json.loads((runs[0] / json_name).read_text())
        _, rows = read_rows(runs[0] / csv_name)
        interval_rows = [r for r in rows if r['value'] == '' and 'correlation' not in r['quantity']]
        assert len(interval_rows) == count_statistics(report) > 0
    manifest = json.loads((runs[0] / 'left_cavity_manifest.json').read_text())
    assert any(str(p).endswith('_cavity_statistics_summary.csv') for p in json.dumps(manifest).split('"'))

    result = compare_regions(runs[0] / 'left_cavity_statistics.json', runs[1] / 'right_cavity_statistics.json',
                             tmp_path / 'comparison')
    _, rows = read_rows(tmp_path / 'comparison' / 'comparison_summary.csv')
    by_quantity = {r['quantity']: r for r in rows}
    assert {'left_volume_A3', 'right_volume_A3', 'paired_volume_difference_A3'} <= set(by_quantity)
    assert float(by_quantity['paired_volume_difference_A3']['mean']) == result['paired_volume_difference_A3']['mean']
    assert len(rows) == count_statistics(result)
