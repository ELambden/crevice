"""Connected transverse free space in an atom-sphere model.

A section is eligible only when its grid component cannot reach the transverse
search boundary. Both nodes and edges exclude probe-inflated atom spheres.
This is a finite-resolution enclosure test, not a continuum proof. A single
axis describes only passages that are monotone along that direction.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from .geometry import add, sub, scale, dot, distance, molecular_basis, principal_axes
from .models import Atom, Coord, ChannelPoint, PoreProfile
from .radii import atom_vdw_radius
from .spatial import SpatialIndex


class ChannelResolutionError(ValueError):
    """A supplied region cannot be resolved as a unique through channel."""


#: The eight in-plane neighbour steps, in the order the flood fill visits them.
_STEPS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


@dataclass
class Section:
    """One enclosed transverse section of the channel at axial coordinate ``t``.

    Attributes
    ----------
    t : float
        Axial coordinate of the plane, Å, relative to the axis origin.
    indices : tuple of tuple of int
        In-plane grid indices ``(i, j)`` of the connected free nodes.
    center : Coord
        Section centre, Å: the node of largest 3D atom-surface clearance (after
        :meth:`SectionGeometry.refine`, the refined local maximum).
    clearance : float
        Atom-surface clearance at ``center``, Å.
    """
    t: float
    indices: tuple[tuple[int, int], ...]
    center: Coord
    clearance: float


class SectionGeometry:
    """Transverse grids perpendicular to a channel axis, and their enclosed sections.

    Each plane at axial coordinate ``t`` carries a square grid of
    ``2*ceil(extent/spacing)+1`` nodes per side spanned by the in-plane basis
    ``(u, v)``. Atoms are van der Waals spheres inflated by ``probe``; a node is
    free outside every inflated disk of that plane.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms.
    origin : Coord
        Axis origin, Å.
    direction : Coord
        Unit axis direction.
    spacing : float, default 0.5
        In-plane grid spacing, Å.
    extent : float, default 6.0
        Half-width of the transverse search square, Å (``search_radius``).
    probe : float, default 0.0
        Enclosure-probe radius, Å, added to every atom radius.
    basis : tuple of Coord, optional
        In-plane axes ``(u, v)``; default from
        :func:`crevice.geometry.molecular_basis`.

    Raises
    ------
    ValueError
        If a plane would need more than 250000 grid points.
    """
    def __init__(self, atoms, origin, direction, *, spacing=0.5, extent=6.0,
                 probe=0.0, basis=None):
        self.atoms = tuple(atoms)
        self.origin, self.direction = origin, direction
        self.u, self.v = basis or molecular_basis(direction, atoms)
        self.spacing, self.probe = spacing, probe
        self.half = int(math.ceil(extent / spacing))
        if (2 * self.half + 1)**2 > 250_000:
            raise ValueError("Section search exceeds 250000 points; reduce search_radius or increase section_spacing")
        self.spatial = SpatialIndex(self.atoms)
        self.projected = [(dot(sub(a.coord, origin), self.u),
                           dot(sub(a.coord, origin), self.v),
                           dot(sub(a.coord, origin), direction), atom_vdw_radius(a)) for a in atoms]

    def position(self, i, j, t):
        """Cartesian position of in-plane node ``(i, j)`` on the plane at ``t``.

        Parameters
        ----------
        i, j : int or float
            In-plane grid indices along ``u`` and ``v`` (0 = on the axis).
        t : float
            Axial coordinate, Å.

        Returns
        -------
        Coord
            ``origin + t*direction + i*spacing*u + j*spacing*v``, Å.
        """
        return add(self.origin, add(scale(self.direction, t),
                   add(scale(self.u, i * self.spacing), scale(self.v, j * self.spacing))))

    def _plane_masks(self, t):
        """Free nodes and blocked grid edges of one plane (numpy, exact).

        Returns ``(free, blocked)``: ``free`` is an ``(n, n)`` boolean array
        (index ``[i+h, j+h]``) of nodes outside every probe-inflated atom disk;
        ``blocked[k]`` marks, for the k-th neighbour step in ``_STEPS``, edges
        from node ``(i, j)`` whose closed segment touches a disk.

        The arithmetic is the same, operation by operation, as the scalar
        tests it replaced (node: ``(i*ds-x)**2 + (j*ds-y)**2 > rr + 1e-12``;
        edge: closest point of the segment to the disk centre), so results are
        bit-identical. Each disk is only evaluated on the index box that can
        reach it (disk radius plus one diagonal step, with a margin); nodes
        and edges outside that box cannot touch the disk. This is a speed-up
        only: it changes no result.
        """
        import numpy as np
        h, ds = self.half, self.spacing
        n = 2*h+1
        grid = np.arange(-h, h+1)*ds
        index = np.arange(-h, h+1)
        inside = np.zeros((n, n), dtype=bool)
        blocked = np.zeros((len(_STEPS), n, n), dtype=bool)
        reach = ds*math.sqrt(2.0)
        for x, y, z, r in self.projected:
            if not abs(z-t) <= r + self.probe:
                continue
            # Exact intersection of the inflated sphere with this plane.
            rr = (r + self.probe)**2 - (z-t)**2
            bound = math.sqrt(max(rr, 0.0) + 1e-12) + 1e-6
            for extent, target in ((bound, "node"), (bound + reach, "edge")):
                i0 = max(-h, int(math.floor((x-extent)/ds)) - 1)
                i1 = min(h, int(math.ceil((x+extent)/ds)) + 1)
                j0 = max(-h, int(math.floor((y-extent)/ds)) - 1)
                j1 = min(h, int(math.ceil((y+extent)/ds)) + 1)
                if i0 > i1 or j0 > j1:
                    continue
                si, sj = slice(i0+h, i1+h+1), slice(j0+h, j1+h+1)
                if target == "node":
                    inside[si, sj] |= ((grid[si]-x)**2)[:, None] + ((grid[sj]-y)**2)[None, :] <= rr + 1e-12
                    continue
                ii = index[si].astype(float)[:, None]
                jj = index[sj].astype(float)[None, :]
                for k, (di, dj) in enumerate(_STEPS):
                    fraction = np.maximum(0.0, np.minimum(1.0, ((x/ds-ii)*di + (y/ds-jj)*dj)/(di*di+dj*dj)))
                    blocked[k, si, sj] |= ((ii+fraction*di)*ds-x)**2 + ((jj+fraction*dj)*ds-y)**2 <= rr + 1e-12
        return ~inside, blocked

    def sections(self, t):
        """Enclosed sections (connected free components) of the plane at ``t``.

        A node is free when the probe-inflated atom disks of this plane leave
        it uncovered; 8-connected free nodes are joined when the closed edge
        between them also stays outside every disk. Components that reach the
        search boundary are exterior and dropped; each remaining component is
        a :class:`Section` centred on its node of largest 3D clearance.
        """
        h = self.half
        free_mask, blocked = self._plane_masks(t)
        free = {(int(i)-h, int(j)-h) for i, j in zip(*free_mask.nonzero())}
        result = []
        while free:
            first = min(free)
            free.remove(first)
            stack, component, exterior = [first], [], False
            while stack:
                node = stack.pop()
                component.append(node)
                i, j = node
                exterior |= abs(i) == h or abs(j) == h
                for k, (di, dj) in enumerate(_STEPS):
                    other = (i+di, j+dj)
                    if other not in free:
                        continue
                    # Minimum distance of each disk center to the closed edge.
                    if not blocked[k, i+h, j+h]:
                        free.remove(other)
                        stack.append(other)
            if not exterior:
                # Real 3D clearance, not the distance to planar disk boundaries.
                scored = [(self.spatial.nearest_surface(self.position(i, j, t))[0], i, j)
                          for i, j in component]
                clearance, i, j = max(scored, key=lambda row: (round(row[0], 10), -row[1]**2-row[2]**2, -row[1], -row[2]))
                result.append(Section(t, tuple(sorted(component)), self.position(i, j, t), clearance))
        return result

    def refine(self, section, steps):
        """Local maximum, reached by collision-free moves inside this component.

        Starting at the section centre, try the eight in-plane moves of half a grid
        spacing, take the one that most increases clearance if the straight move
        clears the probe, halve the step and repeat.

        Parameters
        ----------
        section : Section
            Section to refine.
        steps : int
            Number of halving iterations (``refinement_steps``).

        Returns
        -------
        Section
            Same plane and nodes, with the refined centre and its clearance.
        """
        best, clearance = section.center, section.clearance
        step = self.spacing / 2
        for _ in range(steps):
            improved, score = best, clearance
            for di, dj in ((-1, -1), (-1, 0), (-1, 1), (0, -1),
                           (0, 1), (1, -1), (1, 0), (1, 1)):
                p = add(best, add(scale(self.u, di*step), scale(self.v, dj*step)))
                value = self.spatial.nearest_surface(p)[0]
                if value > score + 1e-12 and self.spatial.segment_clearance(best, p) > self.probe:
                    improved, score = p, value
            best, clearance = improved, score
            step /= 2
        return Section(section.t, section.indices, best, clearance)


def _path_width(geometry, path):
    """Bottleneck of a sampled path: its smallest straight-segment clearance.

    A closed segment contains its end points, so this is also a lower bound
    on every sample clearance. A single-section path has its own clearance.
    """
    if len(path) == 1:
        return path[0].clearance
    return min(geometry.spatial.segment_clearance(a.center, b.center) for a, b in zip(path, path[1:]))


def _runs(geometry, ts, refinement_steps):
    """Longest collision-free paths through adjacent enclosed sections.

    Never interpolate across an absent slice or jump through a wall. A path
    with a sharp bend that cannot be joined is unresolved at this resolution.

    Each enclosed section extends one path ending in the previous plane that
    it can join by a straight segment clearing the enclosure probe. The
    predecessor is ranked by (1) path length, (2) the widest path, i.e. the
    largest bottleneck (minimum segment clearance, including the new segment),
    (3) the shorter lateral step, and (4) section discovery order, which is
    deterministic (smallest transverse grid index first). Rule (2) matters when
    a plane holds, besides the lumen, a small enclosed pocket in the wall (often
    a single grid point where the enclosure probe fits between atoms) that is
    joinable from both neighbouring lumen sections. Ranking by length alone
    tied, and discovery order could then route the profile through the pocket.
    """
    previous, paths = [], []
    for t in ts:
        current = []
        for section in geometry.sections(t):
            section = geometry.refine(section, refinement_steps)
            best, best_key = ([], section.clearance), None
            for path, width in previous:
                edge = geometry.spatial.segment_clearance(path[-1].center, section.center)
                if edge <= geometry.probe:
                    continue
                joined = min(width, edge)
                key = (len(path), round(joined, 10), -round(distance(path[-1].center, section.center), 10))
                if best_key is None or key > best_key:
                    best, best_key = (path, joined), key
            path = best[0] + [section]
            current.append((path, best[1]))
            paths.append(path)
        previous = current
    return paths


def _open_ends(geometry, path):
    """Require unobstructed continuation beyond both ends of the atom envelope."""
    low = min(z-r for _, _, z, r in geometry.projected) - geometry.probe - 1
    high = max(z+r for _, _, z, r in geometry.projected) + geometry.probe + 1
    margins = []
    for section, target in ((path[0], low), (path[-1], high)):
        end = add(section.center, scale(geometry.direction, target-section.t))
        margins.append(geometry.spatial.segment_clearance(section.center, end))
    return margins



def _refine_mouth(geometry, section, sign, step, refinement_steps, tolerance=0.025):
    """Bracket loss of lateral enclosure along a clear exit.

    The bracket describes axial sampling only, excluding transverse-grid and
    enclosure-probe uncertainty. It is a geometric mouth, not a biological label.
    """
    inside = section
    envelope = [z + sign * (r + geometry.probe + 1)
                for _, _, z, r in geometry.projected]
    limit = max(envelope) if sign > 0 else min(envelope)

    def continuation(t):
        choices = [s for s in geometry.sections(t)
                   if geometry.spatial.segment_clearance(inside.center, s.center) > geometry.probe]
        if len(choices) != 1:
            return None
        return geometry.refine(choices[0], refinement_steps)

    outside = inside.t + sign * step
    while sign * (outside-limit) < 0:
        candidate = continuation(outside)
        if candidate is None:
            break
        inside = candidate
        outside += sign * step
    else:
        outside = limit
    while abs(outside-inside.t) > tolerance:
        middle = (outside+inside.t)/2
        candidate = continuation(middle)
        if candidate is None:
            outside = middle
        else:
            inside = candidate
    return inside, {"coordinate_A": (inside.t+outside)/2,
                    "bracket_A": sorted((inside.t, outside)),
                    "enclosed_coordinate_A": inside.t,
                    "open_coordinate_A": outside,
                    "axial_bracket_width_A": abs(outside-inside.t)}


#: Candidate enclosure probes (Å) tried by the automatic choice, in order.
AUTO_ENCLOSURE_LADDER = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20,
                         1.30, 1.40, 1.50, 1.60, 1.80, 2.00, 2.20, 2.40, 2.70, 3.00)
#: Legacy fixed enclosure probe (Å); still used where no probe can be chosen automatically.
LEGACY_ENCLOSURE_RADIUS = 0.8
AUTO_ENCLOSURE_RULE = (
    "smallest probe of the most persistent channel: every ladder probe is resolved with the unchanged "
    "rules; resolved results are grouped into channels (same axis, same exit types, minimum radius within "
    "max(0.05 A, 2 %), axial ranges overlapping >= 60 % of their union); only channels resolved at >= 2 "
    "probes are eligible (a channel seen at one ladder probe is re-tested at the two half-step probes "
    "beside it); the eligible channel resolved at the most probes wins (ties: smaller probe) and its "
    "smallest probe is used; with no eligible channel the profile is unresolved")
_SAME_MIN_RADIUS_A = 0.05
_SAME_MIN_RADIUS_FRACTION = 0.02
_SAME_RANGE_FRACTION = 0.6


def _probe_candidate_record(probe, profile=None, error=None, role="ladder"):
    """One row of ``enclosure_probe_selection.candidates``."""
    if profile is None:
        first = str(error).split(". ")[0].split("; ")[0]
        return {"probe_A": probe, "role": role, "status": "unresolved", "reason": first}
    md = profile.metadata
    return {"probe_A": probe, "role": role, "status": "resolved", "min_radius_A": profile.min_radius,
            "continuous_bottleneck_radius_A": md["continuous_bottleneck_radius_A"],
            "axial_range_A": [profile.points[0].t, profile.points[-1].t],
            "path_type": md.get("path_type", "through_channel")}


def _same_channel(a, b):
    """Whether two resolved profiles (different probes) describe the same channel.

    Same axis direction (|cos| >= 0.99), same exit types per end, minimum
    radii within ``max(0.05 Å, 2 %)`` and axial ranges (sample centres
    projected on ``a``'s axis) overlapping by at least 60 % of their union.
    """
    if abs(dot(a.axis_direction, b.axis_direction)) < 0.99:
        return False
    types = lambda p: tuple(sorted((k, e["type"]) for k, e in (p.metadata.get("exits") or {}).items()))
    if types(a) != types(b):
        return False
    tolerance = max(_SAME_MIN_RADIUS_A, _SAME_MIN_RADIUS_FRACTION*max(a.min_radius, b.min_radius))
    if abs(a.min_radius-b.min_radius) > tolerance:
        return False
    span = lambda p: sorted(dot(sub(q.position, a.axis_origin), a.axis_direction) for q in (p.points[0], p.points[-1]))
    (a0, a1), (b0, b1) = span(a), span(b)
    overlap = max(0.0, min(a1, b1)-max(a0, b0))
    union = max(a1, b1)-min(a0, b0)
    return union > 0 and overlap >= _SAME_RANGE_FRACTION*union


def _auto_enclosure_profile(atoms, *, probe_radius, origin, **kwargs):
    """Choose the enclosure probe automatically and return that profile.

    See ``connected_profile`` (``enclosure_radius=None``) and
    ``AUTO_ENCLOSURE_RULE``. Every candidate is resolved by
    ``connected_profile`` with an explicit probe, so the returned profile is
    exactly what that explicit probe gives, plus
    ``metadata["enclosure_probe_selection"]``.
    """
    ladder = sorted({max(float(probe_radius), q) for q in AUTO_ENCLOSURE_LADDER})
    origin_clearance = SpatialIndex(atoms).nearest_surface(origin)[0] if origin is not None else math.inf
    if origin_clearance <= ladder[0]:
        raise ValueError("Channel origin must have clearance greater than the enclosure and measurement probes "
                         f"(automatic enclosure probe: the smallest candidate is {ladder[0]:.2f} A)")
    rows, results = [], {}

    def attempt(probe, role):
        if probe >= origin_clearance:
            rows.append({"probe_A": probe, "role": role, "status": "not_tried",
                         "reason": "at or above the clearance of the supplied origin"})
            return None
        try:
            profile = connected_profile(atoms, probe_radius=probe_radius, origin=origin,
                                        enclosure_radius=probe, **kwargs)
        except ChannelResolutionError as exc:
            rows.append(_probe_candidate_record(probe, error=exc, role=role))
            return None
        rows.append(_probe_candidate_record(probe, profile, role=role))
        results[probe] = profile
        return profile

    for probe in ladder:
        attempt(probe, "ladder")

    def grouped():
        groups = []
        for probe in sorted(results):
            for group in groups:
                if _same_channel(results[group[0]], results[probe]):
                    group.append(probe)
                    break
            else:
                groups.append([probe])
        return groups

    groups = grouped()
    if groups and max(len(g) for g in groups) < 2:
        # Only isolated ladder probes resolved: re-test half a ladder step on
        # each side before trusting (or refusing) a channel.
        for group in groups:
            k = ladder.index(group[0])
            for other in (ladder[k-1] if k > 0 else None, ladder[k+1] if k+1 < len(ladder) else None):
                if other is not None:
                    attempt(round((group[0]+other)/2, 6), "stability_check")
        groups = grouped()
    group_of = {probe: g for g, members in enumerate(groups, 1) for probe in members}
    for row in rows:
        if row["probe_A"] in group_of:
            row["group"] = group_of[row["probe_A"]]
    record = {"mode": "auto", "rule": AUTO_ENCLOSURE_RULE, "ladder_A": list(AUTO_ENCLOSURE_LADDER),
              "probe_radius_A": probe_radius, "candidates": rows,
              "groups": [{"group": g, "probes_A": members, "support": len(members),
                          "min_radius_A": results[members[0]].min_radius,
                          "continuous_bottleneck_radius_A": results[members[0]].metadata["continuous_bottleneck_radius_A"],
                          "axial_range_A": [results[members[0]].points[0].t, results[members[0]].points[-1].t]}
                         for g, members in enumerate(groups, 1)],
              "criteria": {"same_axis_abs_cos_min": 0.99,
                           "same_min_radius_tolerance": f"max({_SAME_MIN_RADIUS_A} A, {_SAME_MIN_RADIUS_FRACTION:.0%})",
                           "same_axial_range_overlap_of_union_min": _SAME_RANGE_FRACTION,
                           "minimum_support": 2}}
    eligible = [g for g in groups if len(g) >= 2]
    if not eligible:
        tried = [r for r in rows if r["status"] == "unresolved"]
        reasons = {}
        for r in tried:
            reasons.setdefault(r["reason"], []).append(r["probe_A"])
        summary = "; ".join(f"{reason} ({', '.join(f'{p:.3g}' for p in probes)} A)" for reason, probes in reasons.items())
        isolated = ", ".join(f"{g[0]:.3g} A (minimum radius {results[g[0]].min_radius:.3f} A)" for g in groups)
        record.update(chosen_A=None, stability=None,
                      reason=("only isolated probes resolved a channel" if groups
                              else "no candidate probe resolved a channel"))
        message = (f"No enclosure probe resolved a stable channel (automatic choice, "
                   f"{sum(r['status'] != 'not_tried' for r in rows)} probes tried, "
                   f"{ladder[0]:.2f}-{ladder[-1]:.2f} A).")
        if groups:
            message += (f" A channel resolved only at isolated probes, not at a neighbouring probe: {isolated}; "
                        "it is not reported automatically.")
        if summary:
            message += f" Failures: {summary}."
        error = ChannelResolutionError(message + " An explicit enclosure_radius (--enclosure-radius) runs a single probe")
        error.probe_selection = record
        raise error
    best = max(eligible, key=lambda members: (len(members), -members[0]))
    chosen = best[0]
    record.update(chosen_A=chosen, chosen_group=group_of[chosen], stability="plateau",
                  reason=(f"smallest probe of channel group {group_of[chosen]}, which resolved at "
                          f"{len(best)} of {sum(r['status'] != 'not_tried' for r in rows)} probes tried"))
    profile = results[chosen]
    return replace(profile, metadata={**profile.metadata, "enclosure_probe_selection": record})


def connected_profile(atoms, *, axis, origin, samples, search_radius, refinement_steps,
                      probe_radius, section_spacing=0.5, enclosure_radius=0.8,
                      lateral_exits=False, exit_bulk_radius=6.0, exit_spacing=0.5):
    """Resolve one connected, laterally enclosed through-channel and sample it.

    A section is a connected component of transverse grid points where an
    enclosure-probe sphere fits and that cannot reach the search boundary.
    Sections in adjacent planes join only by a straight segment that clears
    the enclosure probe. A plane may hold several sections: the lumen and
    small enclosed pockets in the wall. Pockets are not discarded by size,
    because a narrow lumen can itself be a single grid point; instead the path
    rule chooses the lumen. Each section extends the longest joinable path;
    equal-length alternatives are ranked by the widest path (largest
    bottleneck, i.e. minimum straight-segment clearance), then the shorter
    lateral step, then section discovery order (see ``_runs``). After
    resampling between the refined mouths, full-length paths that start at the
    lower mouth are ranked by bottleneck, then by the distance of their middle
    sample to the origin, then by coordinates. Every choice is deterministic.
    If no path satisfies the rules the channel stays unresolved
    (``ChannelResolutionError``); nothing is interpolated.

    By default both ends must continue straight along the axis, unobstructed,
    beyond the atom envelope (a through-channel). With ``lateral_exits=True``
    an end whose straight exit is blocked (a capped end) is accepted only if a
    lateral exit to bulk exists (:func:`crevice.channel_exits.lateral_exit`).
    The axial profile, its samples and ``min_radius`` are unchanged by this
    option; each end is reported in ``metadata["exits"]`` (``"axial"``, or
    ``"lateral"`` with every distinct leg in ``legs``), ``metadata["path_type"]``
    is ``"through_channel"`` or ``"capped_channel_lateral_exit"``, and
    ``path_bottleneck_radius_A`` is the smaller of the axial continuous
    bottleneck and the widest leg bottleneck of each capped end.

    ``enclosure_radius=None`` chooses the enclosure probe automatically:
    every probe of ``AUTO_ENCLOSURE_LADDER`` (0.5-3.0 Å) is resolved with the
    rules above, results are grouped into channels, and the smallest probe of
    the channel resolved at the most probes (at least two) is used
    (``AUTO_ENCLOSURE_RULE``; design and justification in the ``profile``
    tool guide). The profile is then exactly the one that probe gives
    explicitly, with ``metadata["enclosure_probe_selection"]`` recording the
    ladder, every candidate's outcome, the groups, the chosen probe and the
    reason. If no channel resolves at two probes, ``ChannelResolutionError``
    names the probes and failure reasons and carries the record as
    ``probe_selection``. An explicit number runs that probe only and adds no
    selection metadata.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (already selected).
    axis : str or Coord
        ``"auto"`` (try the principal directions and keep a unique resolved
        channel), ``"x"``, ``"y"``, ``"z"`` or a direction vector.
    origin : Coord or None
        Seed inside the intended channel, Å; it must clear every atom by more
        than both probes. ``None`` uses the atom centroid.
    samples : int
        Minimum number of output samples; more are added so that axial steps
        stay at most 0.75 Å.
    search_radius : float
        Half-width of the transverse search square, Å.
    refinement_steps : int
        Section-centre refinement iterations (:meth:`SectionGeometry.refine`).
    probe_radius : float
        Measurement probe radius, Å, subtracted from each reported radius
        (``raw_clearance`` keeps the unsubtracted value).
    section_spacing : float, default 0.5
        Transverse grid spacing, Å.
    enclosure_radius : float or None, default 0.8
        Enclosure probe, Å: it decides which wall gaps count as closed and which
        sections connect; it does not change the measured radii. ``None`` =
        automatic choice (above). This low-level function defaults to 0.8 Å;
        :func:`crevice.pore_profile` and the command line default to automatic.
    lateral_exits : bool, default False
        Accept capped ends that have a lateral exit (above).
    exit_bulk_radius : float, default 6.0
        With ``lateral_exits``: rolling bulk-probe radius, Å, that defines bulk
        solvent; exits more than twice this apart are distinct legs.
    exit_spacing : float, default 0.5
        With ``lateral_exits``: lattice spacing, Å, of the exit-path search.

    Returns
    -------
    PoreProfile
        Method ``"axial-connected"``: refined, resampled centre-line samples
        (``radius`` = clearance minus ``probe_radius``; ``t`` = axial coordinate,
        Å) with metadata recording the axis, mouths, enclosure probe, section
        basis and search radius, continuous bottleneck, exits and path type.

    Raises
    ------
    ChannelResolutionError
        When no unique connected through-channel (or accepted capped channel)
        is found at this resolution.
    ValueError
        For an origin inside an atom or invalid settings.
    """
    if enclosure_radius is None:
        return _auto_enclosure_profile(
            atoms, axis=axis, origin=origin, samples=samples, search_radius=search_radius,
            refinement_steps=refinement_steps, probe_radius=probe_radius,
            section_spacing=section_spacing, lateral_exits=lateral_exits,
            exit_bulk_radius=exit_bulk_radius, exit_spacing=exit_spacing)
    from .geometry import axis_from_name
    inferred_origin, explicit, label = axis_from_name(axis, atoms)
    supplied_seed = origin is not None
    origin = inferred_origin if origin is None else origin
    if supplied_seed and SpatialIndex(atoms).nearest_surface(origin)[0] <= max(probe_radius, enclosure_radius):
        raise ValueError("Channel origin must have clearance greater than the enclosure and measurement probes")
    automatic = isinstance(axis, str) and axis.lower() == "auto"
    directions = principal_axes(atoms) if automatic else (explicit,)
    candidates, diagnostics = [], []
    for direction in directions:
        geometry = SectionGeometry(atoms, origin, direction, spacing=section_spacing,
                                   extent=search_radius, probe=max(probe_radius, enclosure_radius))
        projections = [z for _, _, z, _ in geometry.projected]
        low, high = min(projections), max(projections)
        # Discovery resolution is independent of requested output point count.
        n = max(3, int(math.ceil((high-low) / 0.75)) + 1)
        ts = [low+(high-low)*i/(n-1) for i in range(n)]
        paths = _runs(geometry, ts, refinement_steps)
        eligible = []
        for path in paths:
            if len(path) < 3:
                continue
            if supplied_seed:
                near = min(path, key=lambda s: abs(s.t))
                if not path[0].t <= 0 <= path[-1].t or geometry.spatial.segment_clearance(origin, near.center) <= geometry.probe:
                    continue
            margins = _open_ends(geometry, path)
            if min(margins) <= geometry.probe and not lateral_exits:
                continue
            # Prioritize axial extent, with clearance as a deterministic tie.
            span = path[-1].t-path[0].t
            eligible.append((span, min(p.clearance for p in path), path, margins))
        best = max(eligible, key=lambda row: row[:2], default=None)
        if best:
            # A short enclosed fragment is not a full channel when another
            # reachable enclosed section lies beyond a gap in its assignment.
            path = best[2]
            beyond = [p[-1] for p in paths if p[-1].t < path[0].t-1e-8 or p[-1].t > path[-1].t+1e-8]
            if lateral_exits:
                # Beyond a capped end, enclosed space (e.g. under the cap) is
                # the region a lateral exit may use, not a channel continuation.
                capped = [margin <= geometry.probe for margin in best[3]]
                beyond = [s for s in beyond if not capped[0 if s.t < path[0].t else 1]]
            if any(geometry.spatial.segment_clearance(
                    path[0 if section.t < path[0].t else -1].center, section.center) > geometry.probe
                   for section in beyond):
                best = None
        diagnostics.append({"direction": direction, "enclosed_sections": len(paths),
                            "through_span_A": best[0] if best else 0.0})
        if best:
            mid = best[2][len(best[2])//2]
            ambiguous = not supplied_seed and any(row[0] >= 0.9*best[0] and
                distance(min(row[2], key=lambda s: abs(s.t-mid.t)).center, mid.center) > 2*section_spacing
                for row in eligible)
            candidates.append((best[0], best[1], geometry, best[2], ambiguous))
    candidates.sort(key=lambda row: row[:2], reverse=True)
    if not candidates:
        raise ChannelResolutionError("No connected channel with two open ends was resolved. Supply a channel axis/origin, "
                         "increase search_radius, reduce section_spacing or adjust enclosure_radius; use cavity analysis for closed voids. "
                         "search_radius=0 requests an unvalidated fixed-axis clearance scan.")
    if len(candidates) > 1 and candidates[1][0] >= 0.9*candidates[0][0]:
        raise ChannelResolutionError("Ambiguous channel axis: multiple directions have comparable connected spans; supply an explicit axis and origin")
    _, _, geometry, coarse, ambiguous = candidates[0]
    if ambiguous:
        raise ChannelResolutionError("Ambiguous channel region: multiple disconnected passages; supply an origin inside the intended channel")
    step = coarse[1].t-coarse[0].t
    lower, lower_mouth = _refine_mouth(geometry, coarse[0], -1, step, refinement_steps)
    upper, upper_mouth = _refine_mouth(geometry, coarse[-1], 1, step, refinement_steps)
    low, high = lower.t, upper.t
    capped = [margin <= geometry.probe for margin in _open_ends(geometry, [lower, upper])] if lateral_exits else [False, False]
    mouth_adjustments = []
    while True:
        # Keep steps <= discovery spacing even when only a few output samples were requested.
        count = max(samples, int(math.ceil((high-low)/0.75))+1)
        ts = [low+(high-low)*i/(count-1) for i in range(count)]
        runs = _runs(geometry, ts, refinement_steps)
        from_lower = [p for p in runs if abs(p[0].t-low) < 1e-9
                      and distance(p[0].center, lower.center) <= 2*section_spacing]
        reached = max((len(p) for p in from_lower), default=0)
        paths = [p for p in from_lower if len(p) == count
                 and (lateral_exits or min(_open_ends(geometry, p)) > geometry.probe)]
        if paths or not any(capped) or len(mouth_adjustments) >= 8:
            break
        # A capped mouth ends where lateral enclosure is first lost. When the
        # denser resampling finds an unenclosed plane inside the bracketed run
        # (a portal that opens and closes along the axis), move that mouth back
        # to the first loss and resample. Never past the seed plane; axial
        # (open) ends are never moved.
        if capped[1] and 1 <= reached < count and ts[reached-1] >= 0:
            inside = max((p for p in from_lower if len(p) == reached), key=lambda p: _path_width(geometry, p))[-1]
            upper, upper_mouth = _refine_mouth(geometry, inside, 1, ts[reached]-inside.t, refinement_steps)
            if upper.t >= high-1e-9:
                break
            mouth_adjustments.append({"end": "upper", "from_A": high, "to_A": upper.t})
            high = upper.t
            continue
        to_upper = [p for p in runs if abs(p[-1].t-high) < 1e-9]
        start = min((count-len(p) for p in to_upper), default=count)
        if capped[0] and 1 <= start < count and ts[start] <= 0:
            inside = max((p for p in to_upper if count-len(p) == start), key=lambda p: _path_width(geometry, p))[0]
            lower, lower_mouth = _refine_mouth(geometry, inside, -1, inside.t-ts[start-1], refinement_steps)
            if lower.t <= low+1e-9:
                break
            mouth_adjustments.append({"end": "lower", "from_A": low, "to_A": lower.t})
            low = lower.t
            continue
        break
    if not paths:
        where = (f" (no enclosed section joinable from the lower mouth at s = {ts[reached]:.2f} A)"
                 if reached < count else "")
        raise ChannelResolutionError("Channel continuity failed at the requested sampling" + where +
                                     "; reduce section_spacing or supply a better axis/origin. Openings in the "
                                     "channel wall wider than the enclosure probe also break lateral enclosure; "
                                     "a larger enclosure_radius (still below the lumen radius) closes them")
    # All remaining paths start at the refined lower mouth and span every
    # sample, so they differ only where a plane holds more than one joinable
    # section. Take the widest (largest bottleneck), then the one whose middle
    # sample is nearest the origin, then the lowest coordinates.
    path = min(paths, key=lambda p: (-round(_path_width(geometry, p), 10),
                                     round(distance(p[len(p)//2].center, origin), 10),
                                     tuple(c for s in p for c in s.center)))
    points = []
    for i, section in enumerate(path):
        raw, nearest = geometry.spatial.nearest_surface(section.center)
        points.append(ChannelPoint(i, section.center, max(0.0, raw-probe_radius), raw, section.t,
                                   nearest.serial, nearest.residue_key.label))
    segments = [geometry.spatial.segment_clearance(a.position, b.position) for a, b in zip(points, points[1:])]
    exit_metadata = {}
    if lateral_exits:
        exit_metadata = _lateral_exit_metadata(geometry, path, segments, probe_radius,
                                               bulk_radius=exit_bulk_radius, spacing=exit_spacing)
        exit_metadata["exit_policy"]["capped_mouth_adjustments"] = mouth_adjustments
    arc = [0.0]
    for a, b in zip(points, points[1:]):
        arc.append(arc[-1]+distance(a.position, b.position))
    return PoreProfile(origin, geometry.direction, tuple(points), probe_radius, method="axial-connected",
                       metadata={"axis": "auto:connected-sections" if automatic else label,
                                 "axis_candidates": diagnostics, "status": "resolved",
                                 "seed_position": path[len(path)//2].center,
                                 "seed_clearance_A": path[len(path)//2].clearance,
                                 "seed_policy": "interior of the selected connected channel",
                                 "volume_estimate_definition": "integral of inscribed sphere areas; not the full non-circular lumen volume",
                                 "section_spacing_A": section_spacing,
                                 "enclosure_radius_A": geometry.probe,
                                 "section_basis": [geometry.u, geometry.v],
                                 "section_search_radius_A": search_radius,
                                 "section_areas_A2": [len(s.indices)*section_spacing**2 for s in path],
                                 "requested_samples": samples, "samples": count,
                                 "search_radius": search_radius, "refinement_steps": refinement_steps,
                                 "reaction_coordinate": "projection onto selected axis; arc lengths recorded separately",
                                 "arc_lengths_A": arc, "segment_clearances_A": segments,
                                 "continuous_bottleneck_radius_A": min(segments)-probe_radius,
                                 "end_clearances_A": _open_ends(geometry, path),
                                 "channel_mouths": {"lower": lower_mouth, "upper": upper_mouth,
                                     "definition": "loss of connected lateral enclosure along the selected axial exit",
                                     "axial_tolerance_A": 0.025,
                                     "limitations": "brackets exclude transverse-grid and probe uncertainty"},
                                 "selection": "connected laterally enclosed sections with two axial exits",
                                 "limitations": "finite grid; principal-direction proposals; monotone unbranched path; not a unique biological assignment",
                                 **exit_metadata})


def _lateral_exit_metadata(geometry, path, segments, probe_radius, *, bulk_radius, spacing):
    """Classify both ends; report every distinct lateral exit of each capped end, or refuse.

    An end whose straight axial exit clears the enclosure probe is
    ``{"type": "axial"}``. A capped end is ``{"type": "lateral"}`` with all
    distinct legs from :func:`crevice.channel_exits.lateral_exit_legs`
    (ranked by width beyond the mouth) under ``legs``; its
    ``bottleneck_radius_A`` is the largest whole-route leg bottleneck, which
    limits the best route out. No leg: the channel is
    refused.
    """
    from .channel_exits import lateral_exit_legs
    if not (math.isfinite(bulk_radius) and bulk_radius > geometry.probe):
        raise ValueError("exit_bulk_radius must be finite and larger than the enclosure probe")
    if not (math.isfinite(spacing) and spacing > 0):
        raise ValueError("exit_spacing must be finite and positive")
    exits, widths = {}, [min(segments)]
    for label, section, sign, margin in zip(("lower", "upper"), (path[0], path[-1]), (-1, 1),
                                            _open_ends(geometry, path)):
        if margin > geometry.probe:
            exits[label] = {"type": "axial", "straight_exit_clearance_A": margin,
                            "definition": "unobstructed straight continuation along the axis beyond the atom envelope"}
            continue
        legs, info = lateral_exit_legs(geometry, section, sign, end=label, spacing=spacing, bulk_radius=bulk_radius)
        if not legs:
            raise ChannelResolutionError(
                f"The {label} end is capped (straight exit clearance {margin:.2f} A) and no lateral exit "
                "to bulk clears the enclosure probe; the channel is not resolved")
        records = [leg.to_dict(probe_radius) for leg in legs]
        widest = max(legs, key=lambda leg: leg.bottleneck_clearance)
        exits[label] = {"type": "lateral", "end": label, "straight_exit_clearance_A": margin,
                        "leg_count": len(records),
                        "bottleneck_radius_A": max(0.0, widest.bottleneck_clearance-probe_radius),
                        "bottleneck_raw_clearance_A": widest.bottleneck_clearance,
                        **info, "legs": records}
        widths.append(widest.bottleneck_clearance)
    lateral = [k for k, e in exits.items() if e["type"] == "lateral"]
    return {"exits": exits,
            "path_type": "capped_channel_lateral_exit" if lateral else "through_channel",
            "path_bottleneck_radius_A": max(0.0, min(widths)-probe_radius),
            "path_bottleneck_definition": ("smallest of the axial continuous bottleneck and, for each capped end, "
                                           "the bottleneck of its widest lateral-exit leg"),
            "exit_policy": {"lateral_exits": True, "bulk_radius_A": bulk_radius, "exit_grid_spacing_A": spacing,
                            "capped_ends": lateral}}


def channel_cast(frame, profile, *, spacing, min_radius, max_grid_points, max_components,
                 max_export_points, min_component_points, focus_points, focus_radius,
                 min_component_volume, max_component_volume, include_hydrogen, include_hetero):
    """Probe-swept lumen from all connected enclosed probe-center sections.

    The full non-circular section contributes, rather than just the inscribed
    sphere at its center. Sweeping the enclosure probe recovers near-wall volume;
    nodes still satisfy the requested atom clearance. Axial ends are explicit
    mouth planes, so the sweep never extends into bulk solvent beyond them.

    Algorithm: on a 3D grid aligned with the profile's axis and section basis,
    each plane between the profile mouths contributes its enclosed section
    nearest the interpolated profile centre (if the straight segment to it
    clears the probe); these core nodes are swept by a sphere of radius
    ``max(enclosure probe, min_radius) - min_radius``; swept nodes with
    atom-surface clearance >= ``min_radius`` are clustered with continuous edge
    tests (:class:`crevice.grid.GridConnectivity`), and the final components are
    selected by size. Planes with no eligible section are recorded as
    unresolved (and warned about), never filled.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile.
    profile : PoreProfile
        Resolved ``"axial-connected"`` profile from :func:`connected_profile`;
        its atom selection must match ``include_hydrogen``/``include_hetero``.
    spacing : float
        Cast grid spacing, Å (never coarsened automatically).
    min_radius : float
        Minimum atom-surface clearance of a cast node, Å.
    max_grid_points : int
        Largest allowed grid; larger requests raise ``ValueError``.
    max_components : int or None
        Maximum number of components kept.
    max_export_points : int
        Largest number of display samples; the measurement keeps every node.
    min_component_points : int
        Components with fewer nodes are dropped.
    focus_points : sequence of Coord or None
        Keep only nodes within ``focus_radius`` of one of these points.
    focus_radius : float
        Focus crop radius, Å.
    min_component_volume, max_component_volume : float or None
        Inclusive volume limits (Å³) of the final components.
    include_hydrogen, include_hetero : bool
        Atom selection; must equal the profile's.

    Returns
    -------
    VoidCast
        Mode ``"channel"``; metadata include the grid origin, basis and
        dimensions, enclosure and sweep radii, ``export_mode
        = "connected_section_fill"``, ``mouth_planes_A``,
        ``unresolved_section_coordinates_A`` and ``status`` (``resolved`` or
        ``partial_or_empty``).

    Raises
    ------
    ValueError
        For an atom-selection mismatch or a grid above ``max_grid_points``.
    """
    from bisect import bisect_right
    from itertools import product
    from .grid import GridConnectivity
    from .models import VoidCast, VoidComponent
    from .voids import _cluster_point_indices, _select_final_components, _renumber_points, _sample_points_by_axis
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    for name, value in (("include_hydrogen", include_hydrogen), ("include_hetero", include_hetero)):
        if profile.metadata.get(name) != value:
            raise ValueError("Channel profile and cast must use the same atom selection")
    q = max(profile.metadata["enclosure_radius_A"], min_radius)
    geometry = SectionGeometry(atoms, profile.axis_origin, profile.axis_direction, spacing=spacing,
                               extent=profile.metadata["section_search_radius_A"], probe=q,
                               basis=profile.metadata["section_basis"])
    low, high = profile.points[0].t, profile.points[-1].t
    nz = int(math.floor((high-low)/spacing + 1e-9))+1
    h = geometry.half
    dims = (2*h+1, 2*h+1, nz)
    count = math.prod(dims)
    if count > max_grid_points:
        raise ValueError(f"Channel grid needs {count} points at {spacing:g} A; increase max_grid_points "
                         "or explicitly choose a coarser spacing. No automatic coarsening was applied.")
    origin = geometry.position(-h, -h, low)
    basis = (geometry.u, geometry.v, geometry.direction)
    connectivity = GridConnectivity(geometry.spatial, origin, spacing, min_radius, basis=basis)
    ts = [p.t for p in profile.points]
    core, gaps = set(), []
    for k in range(nz):
        t = low+k*spacing
        j = max(1, min(len(ts)-1, bisect_right(ts, t)))
        a, b = profile.points[j-1], profile.points[j]
        f = (t-a.t)/(b.t-a.t)
        target = add(a.position, scale(sub(b.position, a.position), f))
        sections = geometry.sections(t)
        eligible = [s for s in sections if geometry.spatial.segment_clearance(target, s.center) >= q]
        if not eligible:
            gaps.append(t)
            continue
        section = min(eligible, key=lambda s: distance(target, s.center))
        core.update((i+h, j+h, k) for i, j in section.indices)
    # Keep only full 3D components that attach to the measured centerline.
    # Separate planes cannot bypass the continuous segment rule.
    radius = max(0.0, q-min_radius)
    cells = int(math.ceil(radius/spacing))
    offsets = [d for d in product(range(-cells, cells+1), repeat=3)
               if spacing**2*sum(x*x for x in d) <= radius**2+1e-12]
    swept = {tuple(i+d for i, d in zip(idx, offset)) for idx in core for offset in offsets}
    swept = {idx for idx in swept if all(0 <= i < n for i, n in zip(idx, dims))}
    nodes = {}
    for idx in sorted(swept):
        position = connectivity.position(idx)
        raw, nearest = geometry.spatial.nearest_surface(position)
        if raw >= min_radius:
            nodes[idx] = ChannelPoint(0, position, raw, raw, low+idx[2]*spacing,
                                      nearest.serial, nearest.residue_key.label)
    if focus_points:
        nodes = {idx: p for idx, p in nodes.items()
                 if any(distance(p.position, f) <= focus_radius for f in focus_points)}
    clusters = _cluster_point_indices(set(nodes), edge_allowed=lambda a, b:
                                      connectivity.allows(a, b, nodes[a].raw_clearance, nodes[b].raw_clearance))
    candidates = []
    for idx, cluster in enumerate(clusters, 1):
        if len(cluster) < min_component_points:
            continue
        points = _renumber_points(sorted((nodes[i] for i in cluster), key=lambda p: (p.t, p.position)))
        center = tuple(sum(p.position[j] for p in points)/len(points) for j in range(3))
        candidates.append(VoidComponent(idx, "section_channel", center, len(points)*spacing**3, points,
                          nearest_residues=tuple(sorted({p.nearest_residue for p in points})),
                          metadata={"source_component_id": idx, "grid_origin": origin,
                                    "grid_basis": basis, "grid_spacing": spacing,
                                    "source_kind": "probe_enclosed_channel",
                                    "selection_cropped": bool(focus_points)}))
    components, selection = _select_final_components(candidates, max_components=max_components,
                                                     min_volume=min_component_volume, max_volume=max_component_volume)
    full = [p for c in components for p in c.points]
    displayed = _renumber_points(_sample_points_by_axis(full, max_export_points))
    if gaps:
        import warnings
        warnings.warn(f"Channel cast is unresolved in {len(gaps)} sections at {spacing:g} A; refine the cast grid. "
                      "See unresolved_section_coordinates_A.", RuntimeWarning, stacklevel=2)
    return VoidCast(components, displayed, spacing, min_radius, "channel", metadata={
        "axis": profile.metadata["axis"], "axis_origin": profile.axis_origin,
        "axis_direction": profile.axis_direction, "grid_origin": origin, "grid_basis": basis,
        "grid_dimensions": dims, "grid_point_count": count,
        "requested_spacing_A": spacing, "effective_spacing_A": spacing, "spacing_adjusted": False,
        "enclosure_radius_A": q, "probe_sweep_radius_A": radius,
        "export_mode": "connected_section_fill", "measurement_point_count": len(full),
        "core_point_count": len(core), "display_subsampled": len(full) > len(displayed),
        "max_export_points": max_export_points, "max_components": max_components,
        "selection": selection, "focus_points": focus_points, "focus_radius": focus_radius if focus_points else None,
        "min_component_volume": min_component_volume, "max_component_volume": max_component_volume,
        "include_hydrogen": include_hydrogen, "include_hetero": include_hetero,
        "connectivity": connectivity.to_dict(), "unresolved_section_coordinates_A": gaps,
        "status": "resolved" if not gaps and len(components) == 1 else "partial_or_empty",
        "mouth_planes_A": (low, low+(nz-1)*spacing),
        "channel_mouths": profile.metadata.get("channel_mouths"),
        "selection_note": "Probe-swept connected lumen, capped at mouth planes; finite-resolution atom-sphere geometry",
        "boundary_validation": "voxel boundaries not continuously validated",
        "limitations": "Narrower crevices than the enclosure probe are excluded; grid refinement is required to assess volume convergence"})
