"""Synthetic command cases for output-inventory characterisation (software checks only).

Each case builds its own synthetic inputs in a scratch directory, runs one
``crevice`` command through :func:`crevice.cli.main` and returns the sorted
relative paths it wrote plus a few key numbers. The same cases are used by
``tests/test_output_policy.py`` and by the inventory characterisation run for
the 0.3 output policy; the 0.2.0 inventories they produced
are frozen in ``tests/data/legacy_inventory_0.2.0.json``. No case is a
biological result.
"""
from __future__ import annotations

import json
from pathlib import Path


from crevice.cli import main
from crevice.presentation import write_viewer_structure

#: Fast profile settings shared by the synthetic channel cases.
PROFILE_ARGS = ['--offline', '--samples', '21', '--search-radius', '8']
#: Small hydration settings (accepted by every command that has them in 0.2.0).
HYDRATION_ARGS = ['--hydration-sasa-points', '32']


def channel_pdb(root: Path) -> Path:
    """Open synthetic rectangular channel (resolves a through-channel profile)."""
    from test_channel_coordinate import rectangular_channel
    path = root / 'channel.pdb'
    if not path.exists():
        write_viewer_structure(rectangular_channel(), path)
    return path


def pocket_pdb(root: Path) -> Path:
    """Capped synthetic channel: no through-channel, a rolling-probe pocket."""
    from test_channel_coordinate import rectangular_channel
    path = root / 'pocket.pdb'
    if not path.exists():
        write_viewer_structure(rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True), path)
    return path


def md_inputs(root: Path):
    """Tiny 3-frame MD fixture with one water and its reference cavity map."""
    from test_hydration import md_fixture, reference_dx
    root.mkdir(parents=True, exist_ok=True)
    source, trajectory = md_fixture(root)
    return source, trajectory, reference_dx(root / 'ref.dx')


def region_definition(root: Path) -> Path:
    """Prepared synthetic region definition (``region-prepare`` already run)."""
    from test_region_definition import definition
    from crevice.region_definition import prepare_region
    root.mkdir(parents=True, exist_ok=True)
    path, _, _ = definition(root)
    return Path(prepare_region(path, root / 'prepared')['definition_json'])


def inventory(root: Path) -> list[str]:
    """Sorted relative paths of every file under ``root``."""
    root = Path(root)
    return sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file())


def _numbers(out: Path, prefix: str) -> dict:
    numbers = {}
    manifest = out / f'{prefix}_manifest.json'
    if manifest.is_file():
        data = json.loads(manifest.read_text(encoding="utf-8"))
        if 'void_cast' in data:
            numbers['void_cast_volume_A3'] = data['void_cast']['total_volume']
        status = data.get('profile_status')
        if isinstance(status, dict):
            numbers['profile_status'] = status.get('status')
    profile = out / f'{prefix}_profile.json'
    if profile.is_file():
        numbers['min_radius_A'] = json.loads(profile.read_text(encoding="utf-8"))['min_radius']
    frames = out / f'{prefix}_cavity_frames.csv'
    if frames.is_file():
        import csv
        with frames.open(newline='', encoding='utf-8') as handle:
            numbers['cavity_volumes_A3'] = [float(r['volume_A3']) for r in csv.DictReader(handle)]
    hydration = out / f'{prefix}_hydration.json'
    if hydration.is_file():
        data = json.loads(hydration.read_text(encoding="utf-8"))
        numbers['hydration_residue_count'] = data.get('residue_count')
    trajectory = out / f'{prefix}_frames.csv'
    if trajectory.is_file():
        import csv
        with trajectory.open(newline='', encoding='utf-8') as handle:
            numbers['frame_min_radius_A'] = [r['min_radius_A'] for r in csv.DictReader(handle)]
    fingerprint = out / f'{prefix}_fingerprint.json'
    if fingerprint.is_file():
        numbers['g_macro_pS'] = json.loads(fingerprint.read_text(encoding="utf-8"))['conductance']['g_macro_pS']
    evidence = out / f'{prefix}_residue_evidence.json'
    if evidence.is_file():
        data = json.loads(evidence.read_text(encoding="utf-8"))
        numbers['lining_residues'] = len(data.get('lining_residues', []))
    return numbers


def run_case(name: str, scratch: Path, extra: list[str] | None = None) -> dict:
    """Run case ``name`` with extra CLI arguments ``extra`` in ``scratch``.

    Parameters
    ----------
    name : str
        One of :data:`CASES`.
    scratch : Path
        Empty working directory; inputs go to ``scratch/inputs`` and outputs to
        ``scratch/out``.
    extra : list of str, optional
        Additional arguments appended to the command (for example ``['--all']``).

    Returns
    -------
    dict
        ``exit`` (status), ``files`` (sorted relative paths under the output
        directory) and ``numbers`` (key measurements read from the outputs).
    """
    scratch = Path(scratch)
    inputs = scratch / 'inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    out = scratch / 'out'
    argv, prefix = CASES[name](inputs, out)
    code = main(argv + list(extra or []))
    files = inventory(out) if out.exists() else []
    return {'exit': code, 'files': files, 'numbers': _numbers(out, prefix) if out.exists() else {}}


def _publish(inputs, out):
    return (['publish', str(channel_pdb(inputs)), *PROFILE_ARGS, '--cast-spacing', '.5', '--dpi', '40',
             *HYDRATION_ARGS, '--out-dir', str(out), '--prefix', 'c'], 'c')


def _publish_pocket(inputs, out):
    return (['publish', str(pocket_pdb(inputs)), *PROFILE_ARGS, '--cast-spacing', '.5', '--dpi', '40',
             '--min-component-volume', '30', *HYDRATION_ARGS, '--out-dir', str(out), '--prefix', 'p'], 'p')


def _cast(inputs, out):
    return (['cast', str(pocket_pdb(inputs)), '--offline', '--min-component-volume', '30', '--dpi', '40',
             *HYDRATION_ARGS, '--out-dir', str(out), '--prefix', 'p'], 'p')


def _analyze(inputs, out):
    return (['analyze', str(channel_pdb(inputs)), *PROFILE_ARGS, '--dpi', '40', *HYDRATION_ARGS,
             '--out-dir', str(out), '--prefix', 'c'], 'c')


def _evidence(inputs, out):
    source, _trajectory, dx = md_inputs(inputs)
    return (['residue-evidence', str(source), '--offline', '--volume-dx', str(dx), '--dpi', '40', *HYDRATION_ARGS,
             '--out-dir', str(out)], 'crevice')


def _hydration_static(inputs, out):
    source, _trajectory, _dx = md_inputs(inputs)
    return (['hydration', str(source), '--offline', '--dpi', '40', *HYDRATION_ARGS, '--out-dir', str(out)], 'crevice')


def _hydration_trajectory(inputs, out):
    source, trajectory, _dx = md_inputs(inputs)
    return (['hydration', str(source), '--trajectory', str(trajectory), '--max-frames', '3', '--dpi', '40',
             '--bootstrap-replicates', '200', *HYDRATION_ARGS, '--out-dir', str(out)], 'crevice')


def _cavity_trajectory(inputs, out):
    source, trajectory, dx = md_inputs(inputs)
    return (['cavity-trajectory', str(source), str(trajectory), '--reference-volume-dx', str(dx),
             '--geometry-mode', 'reference-region', '--max-frames', '3', '--workers', '1', '--dpi', '40',
             '--bootstrap-replicates', '200', *HYDRATION_ARGS, '--out-dir', str(out)], 'crevice')


def _region_trajectory(inputs, out):
    definition = region_definition(inputs)
    return (['region-trajectory', str(definition), '--max-frames', '2', '--workers', '1', '--dpi', '40',
             '--bootstrap-replicates', '200', *HYDRATION_ARGS, '--out-dir', str(out)], 'synthetic')


def _region_prepare(inputs, out):
    from test_region_definition import definition
    path, _, _ = definition(inputs)
    return (['region-prepare', str(path), '--out-dir', str(out)], 'synthetic')


def _trajectory(inputs, out):
    pdb = channel_pdb(inputs)
    return (['trajectory', str(pdb), str(pdb), '-o', str(out / 't.json'), '--samples', '21', '--search-radius', '8',
             '--no-align', '--confidence', '0', '--dpi', '40', *HYDRATION_ARGS], 't')


def _fingerprint(inputs, out):
    return (['fingerprint', str(channel_pdb(inputs)), *PROFILE_ARGS, '--dpi', '40', '--out-dir', str(out),
             '--prefix', 'c'], 'c')


#: Case name -> builder returning (argv, output prefix).
CASES = {
    'publish': _publish,
    'publish-unresolved': _publish_pocket,
    'cast': _cast,
    'analyze': _analyze,
    'residue-evidence': _evidence,
    'hydration-static': _hydration_static,
    'hydration-trajectory': _hydration_trajectory,
    'cavity-trajectory': _cavity_trajectory,
    'region-trajectory': _region_trajectory,
    'region-prepare': _region_prepare,
    'trajectory': _trajectory,
    'fingerprint': _fingerprint,
}

#: Cases whose 0.2.0 file set ``--legacy-outputs`` must reproduce exactly.
LEGACY_CASES = tuple(CASES)

