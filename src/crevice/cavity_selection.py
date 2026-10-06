"""Dominant interior cores and atom-clear, variable-radius sphere filling.

These transparent geometric heuristics rank cavity regions, not biological
channels. Requested counts are upper limits applied after qualification.

Used by :func:`crevice.rolling.rolling_probe_cast` with
``selection_mode="dominant"`` (the default of ``crevice cast``). The idea is
to report the one or few *substantial* buried chambers of a protein rather
than every small crevice: a chamber is found as a *core* of high clearance,
and its volume is then filled with the largest atom-clear spheres centred in
that core.

Definitions
    *Buried centres*: inner-probe centres inside the outer envelope that are
    enclosed in at least ``ceil(26 * enclosure_fraction)`` ray directions.
    *Support*: a connected component of buried centres at the probe radius.
    *Core*: a connected component of buried centres with clearance at least
    the core radius ``max(probe_radius, min(2 Å, 0.35 * peak clearance))``.
    *Fill*: the union of spheres centred on a region's core points with
    radius ``min(atom clearance, depth inside the envelope)``; each sphere
    is therefore wholly atom-free and inside the envelope.
"""
from __future__ import annotations
import math


def _chain_fields(atoms, coords, origin, spacing, shape):
    """Per-voxel single-chain ray enclosure.

    Returns the sorted chain IDs, the best single-chain blocked-ray count per
    voxel, the owning chain (1-based, the first chain on ties) and the number
    of chains that block at least 6 rays.
    """
    import numpy as np
    from .rolling import _ray_enclosure
    from .radii import atom_vdw_radius
    chains = sorted({a.chain_id for a in atoms})
    best = np.zeros(shape, dtype=np.uint8)
    owner = np.zeros(shape, dtype=np.int16)
    contributors = np.zeros(shape, dtype=np.uint8)
    for ci, chain in enumerate(chains, 1):
        blocked = np.zeros(shape, dtype=bool)
        for atom, coord in zip(atoms, coords):
            if atom.chain_id != chain:
                continue
            radius = atom_vdw_radius(atom)
            lo = np.maximum(np.ceil((coord-radius-origin)/spacing).astype(int), 0)
            hi = np.minimum(np.floor((coord+radius-origin)/spacing).astype(int)+1, shape)
            if np.any(hi <= lo):
                continue
            box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
            ijk = np.indices(tuple(hi-lo)).transpose(1,2,3,0)+lo
            blocked[box] |= np.sum((origin+ijk*spacing-coord)**2, axis=-1) <= radius**2
        rays = _ray_enclosure(blocked)
        contributors += rays >= 6
        take = rays > best
        owner[take] = ci
        best = np.maximum(best, rays)
    return chains, best, owner, contributors


def _paint_balls(indices, radii, spacing, score, labels, label):
    """Paint the exact sampled union of individually atom-clear spheres.

    A voxel inside several spheres goes to the sphere with the largest margin
    ``radius - distance``; ``labels`` and ``score`` are updated in place.
    """
    import numpy as np
    shape = np.asarray(score.shape)
    for index, radius in zip(indices, radii):
        extent = int(math.ceil(radius/spacing))
        lo, hi = np.maximum(index-extent, 0), np.minimum(index+extent+1, shape)
        box = tuple(slice(int(a), int(b)) for a,b in zip(lo,hi))
        axes = np.ogrid[tuple(slice(int(a-i), int(b-i)) for a,b,i in zip(lo,hi,index))]
        distance = np.sqrt(sum(a*a for a in axes))*spacing
        margin = radius-distance
        take = (margin >= -1e-10) & (margin > score[box]+1e-10)
        score[box][take] = margin[take]
        labels[box][take] = label


def dominant_sphere_fill(*, atoms, coords, clearance, centers, depth, origin,
                         spacing, probe_radius, enclosure_fraction, queries,
                         min_component_volume, core_fraction=.08):
    """Qualify dominant buried cores and fill them with atom-clear spheres.

    Algorithm:

    1. Keep probe centres enclosed in at least ``ceil(26 *
       enclosure_fraction)`` ray directions (atoms inflated by the probe).
    2. Label supports (probe-radius connectivity) and cores (connectivity at
       the core radius; see the module description).
    3. A core qualifies if it has at least ``max(8, ceil(0.05 *
       min_component_volume / spacing**3), ceil(core_fraction * largest core))``
       points; weaker cores are reported as ``weak_core``.
    4. With three or more chains, each core is assigned to a chain when that
       chain alone blocks at least 20 rays at 25 % or more of the core's
       voxels and wins at least 60 % of those votes. Qualified cores are
       grouped by (support, chain). Without chain assignment, several cores in
       one support (lobes joined by a neck) become one region spanning the
       whole support.
    5. With three or more chains, supports not already represented can add a
       ``shared_interface_core`` region when enough of their buried space is
       blocked by at least three chains, none dominating (at most 13 rays);
       see ``interface_rule`` in the report.
    6. Candidates are ranked by core size (ties by support id, then chain) and
       painted in that order; overlaps go to the sphere with the larger margin.
       Voxels with negative clearance or non-positive depth are cleared.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (for chain IDs and radii).
    coords : ndarray, shape (n, 3)
        Atom centres in the grid frame (Å).
    clearance : ndarray
        Atom-surface clearance at every grid node (Å).
    centers : ndarray of bool
        Admissible inner-probe centres.
    depth : ndarray
        Distance (Å) inside the outer envelope, already reduced by
        ``min_depth``.
    origin : ndarray, shape (3,)
        Grid-frame position of node ``(0, 0, 0)`` (Å).
    spacing : float
        Grid spacing (Å).
    probe_radius : float
        Inner probe radius (Å).
    enclosure_fraction : float
        Required fraction of the 26 ray directions.
    queries : _SphereQueries
        Atom-sphere queries in the grid frame.
    min_component_volume : float
        Minimum region volume (Å³); also scales the core-size minimum.
    core_fraction : float, default 0.08
        Minimum core size relative to the largest core, in ``(0, 1]``.

    Returns
    -------
    labels : ndarray of int32
        Region label per grid node (0 = none), numbered by rank from 1.
    center_counts : ndarray
        Number of core points per label (index 0 unused).
    report : dict
        Core radius, qualification threshold, each core's decision, each
        region's support/chain/core ids and the fill, interface and
        interpretation definitions.

    Raises
    ------
    ValueError
        If ``core_fraction`` is not in ``(0, 1]``.
    """
    import numpy as np
    from scipy import ndimage
    from .rolling import _connected, _ray_enclosure
    if not math.isfinite(core_fraction) or not 0 < core_fraction <= 1:
        raise ValueError('core_fraction must be in (0, 1]')
    required_rays = math.ceil(26*enclosure_fraction)
    buried = centers & (_ray_enclosure(clearance < probe_radius) >= required_rays)
    empty = np.zeros(clearance.shape, dtype=np.int32)
    if not buried.any():
        return empty, np.zeros(1, dtype=int), {
            'algorithm':'dominant_core_variable_sphere_fill', 'candidates': [],
            'regions': [], 'qualified_count': 0, 'core_radius_A':None,
            'core_volume_fraction':core_fraction,
            'fill_definition':'union of atom-clear spheres centered on qualified cores; each radius=min(atom clearance, outer-sweep distance)',
            'reason':'no probe centers meet the envelope and core-enclosure criteria',
        }
    support = _connected(buried, clearance, origin, spacing, probe_radius, queries)
    peak = float(clearance[buried].max())
    core_radius = max(probe_radius, min(2., .35*peak))
    core = _connected(buried & (clearance >= core_radius), clearance,
                      origin, spacing, core_radius, queries)
    counts = np.bincount(core.ravel())
    minimum = max(8, int(math.ceil(min_component_volume*.05/spacing**3)),
                  int(math.ceil(core_fraction*counts[1:].max())))
    qualified = [i for i in range(1,len(counts)) if counts[i] >= minimum]
    boxes = ndimage.find_objects(core)
    chains = sorted({a.chain_id for a in atoms})
    fields = _chain_fields(atoms, coords, origin, spacing, clearance.shape) if len(chains) >= 3 else None
    groups, decisions = {}, []
    for label in range(1,len(counts)):
        if not counts[label]:
            continue
        decisions.append({'core_id': label, 'core_volume_A3': int(counts[label])*spacing**3,
                          'decision': 'qualified' if label in qualified else 'weak_core'})
        if label not in qualified:
            continue
        box = boxes[label-1]
        idx = np.argwhere(core[box] == label)+np.array([s.start for s in box])
        # A higher-clearance core lies in one small-probe support component.
        parent = int(support[tuple(idx[0])])
        chain = 0
        if fields:
            _, best, owner, _ = fields
            votes = owner[tuple(idx.T)].copy()
            votes[best[tuple(idx.T)] < 20] = 0
            tally = np.bincount(votes, minlength=len(chains)+1)
            winner = int(np.argmax(tally[1:]))+1
            if tally[winner] >= .25*len(idx) and tally[winner] >= .6*tally[1:].sum():
                chain = winner
        # Do not split a two-chain, end-to-end pore such as gramicidin.
        # Multiple strong subunit interiors in larger multimers may remain
        # separate regions even when their mouths share a support component.
        groups.setdefault((parent,chain), []).append((label,idx))
    candidates = []
    represented = set()
    for (parent,chain), rows in groups.items():
        represented.add(parent)
        idx = np.concatenate([row[1] for row in rows])
        if len(rows) > 1 and not fields:
            # Restore the small-probe support joining lobes across a narrow
            # neck; a neck does not create extra requested cavities.
            idx = np.argwhere(support == parent)
        candidates.append({'indices':idx, 'kind':'dominant_core', 'support_id':parent,
                           'chain_id':chains[chain-1] if chain else None,
                           'core_ids':[row[0] for row in rows], 'core_radius_A':core_radius})
    if fields:
        _, best, _, contributors = fields
        balanced = buried & (contributors >= 3) & (best <= 13)
        # Shared space already belonging to a selected chamber is not another
        # cavity. Only an independent small-probe support is eligible here.
        shared_counts = np.bincount(support[balanced], minlength=int(support.max())+1)
        eligible = [(i,int(n)) for i,n in enumerate(shared_counts)
                    if i and i not in represented and n*spacing**3 >= max(8.,min_component_volume*.25)]
        if eligible:
            strongest = max(n for _,n in eligible)
            for parent,n in eligible:
                if n < .5*strongest:
                    continue
                mask = support == parent
                radius = max(probe_radius, min(2., .35*float(clearance[mask].max())))
                cc = _connected(mask & (clearance >= radius), clearance, origin, spacing, radius, queries)
                sizes = np.bincount(cc.ravel()); label = int(np.argmax(sizes[1:]))+1
                idx = np.argwhere(cc == label)
                candidates.append({'indices':idx, 'kind':'shared_interface_core', 'support_id':parent,
                                   'chain_id':None, 'core_ids':[], 'core_radius_A':radius,
                                   'balanced_support_volume_A3':n*spacing**3})
    # Rank before painting so overlap ownership is independent of a requested
    # output count. Equal scores have deterministic core/support tie-breaks.
    candidates.sort(key=lambda x:(-len(x['indices']),x['support_id'],x['chain_id'] or ''))
    score = np.full(clearance.shape, -np.inf, dtype=np.float32)
    labels = empty
    center_counts = [0]
    rows = []
    for label,candidate in enumerate(candidates, 1):
        idx = candidate.pop('indices')
        # Both functions are 1-Lipschitz. A ball no larger than atom clearance
        # and distance to the outer sweep is wholly inside admissible space.
        radii = np.minimum(clearance[tuple(idx.T)], depth[tuple(idx.T)])
        _paint_balls(idx,radii,spacing,score,labels,label)
        center_counts.append(len(idx))
        rows.append({'region_id':label, **candidate, 'core_points':len(idx),
                     'core_volume_A3':len(idx)*spacing**3})
    labels[(clearance < 0) | (depth <= 0)] = 0
    return labels, np.asarray(center_counts), {
        'algorithm':'dominant_core_variable_sphere_fill',
        'core_radius_A':core_radius, 'core_volume_fraction':core_fraction,
        'minimum_core_points':minimum, 'qualified_count':len(rows),
        'candidates':decisions, 'regions':rows,
        'interface_rule':'independent support; at least 3 chains block >=6 rays each; no chain >13; support >=max(8 A^3, min_volume/4) and >=50% of strongest shared support',
        'fill_definition':'union of atom-clear spheres centered on qualified cores; each radius=min(atom clearance, outer-sweep distance)',
        'interpretation':'geometric cavity-region heuristic; not a biological channel classifier',
    }
