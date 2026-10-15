"""Typed residue–residue interactions without the water analyses.

Since CREVICE 0.3 the hydration analysis (water contacts, SASA, water
bridges, water density) is opt-in, while the chemically typed residue–residue
interactions it used to carry (hydrogen bonds, salt-bridge candidates,
aromatic contacts) and the residue conformations still run by default. This
module computes that residue part on its own.

The residue focus, the chemical topology (built with the same water oxygens,
so hydrogen attachment is unchanged) and the per-frame measurement are
exactly those of the hydration workflow
(:func:`crevice.hydration_extensions.attach_static_extensions` and
:func:`crevice.hydration_trajectory.analyze_hydration_trajectory`); only the
water-specific terms are removed from the result (:func:`residue_part`), so
every residue–residue edge and conformation number equals the one written by
the hydration workflow. Outputs are the ``PREFIX_interactions.json``,
``PREFIX_interaction_*`` and ``PREFIX_residue_conformation*`` files of
:func:`crevice.interaction_analysis.write_interaction_bundle`.
"""
from __future__ import annotations

import json
from pathlib import Path

#: Interaction kinds that involve water; reported only with the hydration analysis.
WATER_EDGE_KINDS = frozenset({"water_bridge_hbond"})

WATER_TERMS_NOTE = ("not measured: water bridges and residue-water hydrogen bonds belong to the hydration "
                    "analysis (--outputs hydration or --all)")


def residue_part(result: dict) -> dict:
    """Remove the water-specific terms from an interaction result.

    Parameters
    ----------
    result : dict
        :meth:`crevice.interaction_analysis.InteractionAnalysis.result` output.

    Returns
    -------
    dict
        A copy without water-bridge edges; residue-water hydrogen-bond counts
        are NaN (unmeasured) and their contacts empty; ``metadata`` records
        ``water_terms``.
    """
    import numpy as np
    out = dict(result)
    out["edges"] = {key: value for key, value in result["edges"].items() if key[0] not in WATER_EDGE_KINDS}
    out["hydrogen_bonded_water_count"] = np.full(np.shape(result["hydrogen_bonded_water_count"]), np.nan)
    out["hydrogen_bonded_water_contacts"] = np.zeros((0, 3), dtype=int)
    out["metadata"] = {**result["metadata"], "water_terms": WATER_TERMS_NOTE}
    return out


def static_residue_interactions(frame, *, source_path=None, volume_path=None, residue_ids=None):
    """Typed residue–residue interactions of one structure (no files written).

    The residue focus follows :func:`crevice.hydration_workflow.static_hydration_bundle`:
    ``residue_ids``, else the boundary residues and nonlocal partners of
    ``volume_path``, else every protein residue; HETATM and water residues are
    not typed. A ``.gro`` source is reloaded with its solvent and box.

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    source_path : str or Path, optional
        Source file (a ``.gro`` triggers reloading the full system).
    volume_path : str or Path, optional
        Measured binary cast DX defining the focus.
    residue_ids : sequence of str, optional
        Explicit focus.

    Returns
    -------
    dict or None
        :func:`residue_part` of the interaction result, or ``None`` when the
        focus is empty.

    Raises
    ------
    ValueError
        For focus residues absent from the structure.
    """
    import numpy as np
    from .hydration import WATER_RESNAMES
    from .interaction_analysis import InteractionAnalysis, chemical_topology_from_frame
    box = None
    if volume_path and residue_ids is None:
        from .residue_evidence import boundary_residue_evidence, read_binary_dx
        grid, origin, deltas = read_binary_dx(volume_path)
        evidence = boundary_residue_evidence(frame, grid, origin, deltas)
        residue_ids = [r['residue'] for r in evidence['residues'] if r['role'] != 'other']
    if source_path and Path(source_path).suffix.lower() == '.gro':
        from .parser import load_structure_report
        from .trajectory import _open_universe
        full, _ = load_structure_report(source_path, model_index=frame.model_index, md_selection='all')
        universe = _open_universe(source_path)
        universe.trajectory[frame.model_index]
        box = universe.dimensions
        frame = full
    known = {a.residue_key.label for a in frame.atoms
             if a.element.upper() not in {'H', 'D'} and a.resname.upper() not in WATER_RESNAMES and not a.hetero}
    if residue_ids is None:
        focus = known
    else:
        missing = set(residue_ids) - {a.residue_key.label for a in frame.atoms}
        if missing:
            raise ValueError('Interaction focus residue identities are not present in the structure: ' + str(sorted(missing)))
        focus = set(residue_ids) & known
    if not focus:
        return None
    waters = [i for i, a in enumerate(frame.atoms) if a.resname.upper() in WATER_RESNAMES and a.element.upper() == 'O'
              and (a.occupancy is None or a.occupancy > 0)]
    chemistry = chemical_topology_from_frame(frame, waters, box=box)
    analysis = InteractionAnalysis(chemistry, focus, 1)
    analysis.add(np.asarray([a.coord for a in frame.atoms]), box=box)
    return residue_part(analysis.result())


def trajectory_residue_interactions(topology, trajectory, *, residue_ids=None, context_json=None,
                                    selection='protein and not name H*', water_selection=None, start=0, stop=None,
                                    stride=1, max_frames=500, pbc='check', progress_callback=None):
    """Typed residue–residue interactions over an MD trajectory (no files written).

    Frames, focus and chemistry follow
    :func:`crevice.hydration_trajectory.analyze_hydration_trajectory`; with
    ``context_json`` (a cavity-trajectory statistics file) the focus is every
    residue that was ever lining or a nonlocal partner and the frames are those
    of the geometry run. No SASA, water contact or density is measured.

    Parameters
    ----------
    topology, trajectory : str or Path
        MD files readable by MDAnalysis.
    residue_ids : sequence of str, optional
        Explicit focus (default: the context residues, else every selected residue).
    context_json : str or Path, optional
        ``PREFIX_cavity_statistics.json`` of a cavity or region trajectory run.
    selection : str, default "protein and not name H*"
        MDAnalysis selection of the protein atoms.
    water_selection : str, optional
        Water selection (default :data:`crevice.hydration.DEFAULT_WATER_SELECTION`).
    start, stop, stride : int
        Frame range (source indices).
    max_frames : int, default 500
        Frame guard.
    pbc : {"check", "none", "unwrap"}, default "check"
        Periodic handling.
    progress_callback : callable, optional
        Called as ``(row, n)`` after each frame.

    Returns
    -------
    result : dict or None
        :func:`residue_part` of the interaction result (``None`` for an empty focus).
    times : numpy.ndarray
        Frame times (ps).
    """
    import numpy as np
    from .hydration import DEFAULT_WATER_SELECTION
    from .interaction_analysis import InteractionAnalysis, chemical_topology_from_universe
    from .trajectory import _element_names, _open_universe, load_trajectory_report
    if context_json is not None:
        path = Path(context_json).resolve(strict=True)
        context = json.loads(path.read_text())
        stem = path.name.removesuffix('_cavity_statistics.json')
        reader = json.loads((path.parent / (stem + '_reader.json')).read_text())
        indices = list(range(*reader['frame_indices']))
        chosen = [f for f in indices if f >= start and (stop is None or f < stop) and (f - start) % stride == 0]
        if not chosen:
            raise ValueError('No context frames match the requested range')
        start, stop = chosen[0], chosen[-1] + 1
        stride = chosen[1] - chosen[0] if len(chosen) > 1 else 1
        selection = reader['selection']
        if residue_ids is None:
            residue_ids = [r['residue'] for r in context['residues'] if r['lining_frames'] or r['nonlining_partner_frames']]
    protein, _ = load_trajectory_report(topology, trajectory, selection=selection, start=start, stop=stop,
                                        stride=stride, max_frames=max_frames, pbc=pbc)
    times = np.asarray(protein.times(), dtype=float)
    first = protein.frames[0]
    known = {a.residue_key.label for a in first.atoms}
    focus = known if residue_ids is None else set(residue_ids)
    if focus - known:
        raise ValueError('Interaction focus must contain valid selected protein residue identities')
    if not focus:
        return None, times
    universe = _open_universe(topology, trajectory)
    atoms = universe.atoms
    selected = universe.select_atoms(selection)
    elements = np.asarray(_element_names(atoms))
    water_group = universe.select_atoms(water_selection or DEFAULT_WATER_SELECTION)
    oxygen_indices = water_group.indices[elements[water_group.indices] == 'O']
    # Same residue ids and order as the hydration model (selected protein atoms of the focus).
    targets = sorted({a.residue_key.label for a in first.atoms if a.residue_key.label in focus})
    universe.trajectory[first.frame_index]
    chemistry = chemical_topology_from_universe(universe, selected, first, oxygen_indices)
    analysis = InteractionAnalysis(chemistry, targets, len(protein.frames))
    for row, frame in enumerate(protein.frames):
        universe.trajectory[frame.frame_index]
        box = None if pbc == 'none' else universe.dimensions
        if pbc != 'none' and box is None:
            raise ValueError('Periodic interactions need cell dimensions; explicitly use pbc=none for nonperiodic data')
        analysis.add(atoms.positions, box=box)
        if progress_callback:
            progress_callback(row, len(protein.frames))
    return residue_part(analysis.result()), times


def write_residue_interactions(result, output_dir, prefix, *, times=None, confidence=.95, block_length=None,
                               replicates=2000, seed=20260914) -> dict:
    """Write the residue-interaction files of :func:`static_residue_interactions` or :func:`trajectory_residue_interactions`.

    Parameters
    ----------
    result : dict or None
        The residue interaction result; ``None`` writes
        ``PREFIX_interactions.json`` with an ``unavailable`` status.
    output_dir : str or Path
        Output directory (created).
    prefix : str
        File-name prefix.
    times : array_like, optional
        Frame times (ps) of a trajectory; ``None`` for one structure.
    confidence : float, default 0.95
        Occupancy interval level (0 = descriptive only).
    block_length : int, optional
        Bootstrap block length in frames.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260914
        Bootstrap seed.

    Returns
    -------
    dict of str to str
        Output keys and paths (see
        :func:`crevice.interaction_analysis.write_interaction_bundle`).
    """
    from .interaction_analysis import write_interaction_bundle
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    if result is None:
        path = root / (prefix + '_interactions.json')
        path.write_text(json.dumps({'status': 'unavailable_no_identified_residues', 'edges': [], 'residues': [],
                                    'interpretation': 'No residue focus was available; not evidence of absent interactions'},
                                   indent=2) + '\n')
        return {'interaction_json': str(path)}
    return write_interaction_bundle(result, root, prefix, times=times, confidence=confidence,
                                    block_length=block_length, replicates=replicates, seed=seed)
