"""Widest-path access tunnels on a clearance grid.

A tunnel here is a path on a regular grid from a start point (for example a
buried active site) to the edge of the padded bounding box, along which a
probe of radius ``min_radius`` can move without touching any atom van der
Waals sphere. Among all such paths, CREVICE keeps the *widest* ones: those
that maximise the bottleneck (the smallest clearance along the path), with
shorter length breaking ties.

Main entry points: :func:`find_tunnels` (command ``crevice tunnels``),
:func:`rank_tunnels`, :func:`cluster_tunnels`, :func:`tunnel_profile`.

Each tunnel point's ``t`` is the cumulative path length in Å from the start,
as for explicit-path profiles; its step number along the path is
``ChannelPoint.index``. (:func:`crevice.analysis.find_tunnels` is the
underlying grid search; it returns the same Å path lengths, and every command
and the public ``crevice.find_tunnels`` use this module's :func:`find_tunnels`.)

Assumptions and limitations
    The search is on a grid of the chosen spacing (2 Å by default), so a
    narrow tunnel between grid points can be missed; the absence of a tunnel
    is not proof that none exists. Grid edges are tested analytically along
    the whole segment, so a reported path does clear every atom sphere. The
    box edge stands in for bulk solvent; membranes, crystal contacts and
    missing residues are not modelled. Tunnels are geometric candidates, not
    demonstrated ligand or solvent pathways.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Sequence

from .analysis import find_tunnels as _find_tunnels
from .geometry import distance, normalize, sub
from .models import PoreProfile, StructureFrame, Tunnel
from .radii import RadiusSet, radii_option

#: Meaning of ``ChannelPoint.t`` on tunnel points, recorded in tunnel profiles.
TUNNEL_T_DEFINITION = "cumulative path length in Å from the tunnel start (index = step number)"


@radii_option
def find_tunnels(frame: StructureFrame, *, radii: RadiusSet | str | None = None, **kwargs) -> tuple[Tunnel, ...]:
    """Find the widest clearance-limited grid paths from a start to the box edge.

    Runs the grid search :func:`crevice.analysis.find_tunnels` (keyword
    arguments are passed through unchanged), whose points already carry the
    cumulative path length in Å as ``t`` (:func:`with_path_length`, applied
    again here; it is idempotent). Nothing else about the tunnels changes.

    Algorithm:

    1. Select atoms (heavy, non-HETATM by default) and build a grid over their
       bounding box padded by ``padding``. If the grid would exceed
       ``max_grid_points``, the spacing is enlarged uniformly (recorded in each
       tunnel's metadata as ``spacing_adjusted``).
    2. Keep grid nodes whose atom-surface clearance is at least ``min_radius``.
    3. Attach the start to the nearest kept node within one cell diagonal
       (``sqrt(3) * spacing``) whose connecting segment also clears
       ``min_radius``.
    4. Run a widest-path (maximin) search over 26-neighbour edges whose whole
       segment clears ``min_radius``; equal bottlenecks prefer the shorter path.
    5. Rank reachable nodes on the box boundary by bottleneck (descending), then
       length (ascending); accept them in turn, skipping ends within
       ``3 * spacing`` of an already accepted end, until ``max_tunnels``.

    Parameters
    ----------
    frame : StructureFrame
    **kwargs
        Passed to :func:`crevice.analysis.find_tunnels`:

        start : Coord, optional
            Start point in Å. Defaults to the centroid of the selected atoms.
        spacing : float, default 2.0
            Requested grid spacing in Å.
        min_radius : float, default 0.8
            Required clearance (probe radius) in Å for nodes and edges. It is
            a threshold only; reported radii are not reduced by it.
        max_tunnels : int, default 3
            Maximum number of tunnels returned.
        padding : float, default 3.0
            Å added around the atom bounding box. The box faces act as bulk.
        max_grid_points : int, default 120000
            Grid-size limit that triggers automatic coarsening.
        include_hydrogen : bool, default False
            Treat hydrogen atoms as obstacles.
        include_hetero : bool, default False
            Treat HETATM atoms as obstacles. Note the default differs from
            :func:`crevice.channels.pore_profile`.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    tuple of Tunnel
        Widest first. Empty if the default start (the centroid) is inside an
        atom or too tight, or if no path reaches the boundary. Each tunnel's
        ``points`` begin with the start itself; ``t`` is the cumulative path
        length in Å from the start (so the last point's ``t`` equals
        ``length``) and ``index`` is the step number.

    Raises
    ------
    ValueError
        If an explicit ``start`` has clearance below ``min_radius``, or for
        invalid ``spacing``/``min_radius`` values.

    See Also
    --------
    rank_tunnels, cluster_tunnels, tunnel_profile

    Examples
    --------
    Paths out of a synthetic cylinder (radius 6 Å, open at both ends):

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.tunnels import find_tunnels
    >>> atoms = [Atom(12 * k + j + 1, "C", "ALA", "A", k + 1,
    ...               6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6),
    ...               -10.5 + 1.5 * k, "C")
    ...          for k in range(15) for j in range(12)]
    >>> tunnels = find_tunnels(StructureFrame(tuple(atoms)), start=(0.0, 0.0, 0.0))
    >>> [(round(t.bottleneck_radius, 2), round(t.length, 1)) for t in tunnels[:2]]
    [(2.95, 15.0), (2.95, 17.0)]
    """

    return tuple(with_path_length(tunnel) for tunnel in _find_tunnels(frame, **kwargs))


def with_path_length(tunnel: Tunnel) -> Tunnel:
    """Return the tunnel with each point's ``t`` set to its path length in Å.

    ``t`` becomes the sum of straight segment lengths from the first point
    (0 at the start). Positions, radii, indices (the step numbers) and every
    other field are unchanged.

    Parameters
    ----------
    tunnel : Tunnel

    Returns
    -------
    Tunnel
    """
    points, travelled = [], 0.0
    for index, point in enumerate(tunnel.points):
        if index:
            travelled += distance(tunnel.points[index - 1].position, point.position)
        points.append(replace(point, t=travelled))
    return replace(tunnel, points=tuple(points))


def rank_tunnels(tunnels: Sequence[Tunnel], *, by: str = "bottleneck_then_length") -> tuple[Tunnel, ...]:
    """Sort tunnels by a named, deterministic policy.

    Parameters
    ----------
    tunnels : sequence of Tunnel
    by : {"bottleneck_then_length", "score", "length"}, default "bottleneck_then_length"
        * ``"bottleneck_then_length"``: widest bottleneck first, then shortest.
        * ``"score"``: highest :attr:`~crevice.models.Tunnel.score` first, then widest, then
          shortest.
        * ``"length"``: shortest first, then widest.

    Returns
    -------
    tuple of Tunnel
        Sorted copy; the sort is stable, so full ties keep input order.

    Raises
    ------
    ValueError
        For an unknown policy.
    """

    if by == "score":
        key = lambda tunnel: (-tunnel.score, -tunnel.bottleneck_radius, tunnel.length)
    elif by == "length":
        key = lambda tunnel: (tunnel.length, -tunnel.bottleneck_radius)
    elif by == "bottleneck_then_length":
        key = lambda tunnel: (-tunnel.bottleneck_radius, tunnel.length)
    else:
        raise ValueError("by must be score, length, or bottleneck_then_length")
    return tuple(sorted(tunnels, key=key))


def tunnel_profile(tunnel: Tunnel) -> PoreProfile:
    """View a tunnel as a :class:`~crevice.models.PoreProfile` for code that expects profiles.

    Parameters
    ----------
    tunnel : Tunnel

    Returns
    -------
    PoreProfile
        ``method="tunnel-path"``, ``axis_origin`` at the tunnel start,
        ``axis_direction`` the unit vector from start to end (``(0, 0, 1)`` for
        a zero-length tunnel) and ``probe_radius=0``, so point radii are the raw
        clearances.

    Notes
    -----
    The points keep the tunnel's ``t`` values: cumulative path length in Å
    from the start for tunnels from :func:`find_tunnels`, recorded as
    ``metadata["t_definition"]``. They are distances along this tunnel's own
    path, not positions on a shared axis, so axial statistics such as
    :func:`crevice.ensemble.profile_distribution` do not pool them.
    """

    direction = normalize(sub(tunnel.end, tunnel.start)) if tunnel.length else (0.0, 0.0, 1.0)
    return PoreProfile(
        axis_origin=tunnel.start,
        axis_direction=direction,
        points=tunnel.points,
        method="tunnel-path",
        metadata={"tunnel_id": tunnel.id, "score": tunnel.score, "t_definition": TUNNEL_T_DEFINITION},
    )


def cluster_tunnels(tunnels: Sequence[Tunnel], *, endpoint_cutoff: float = 5.0) -> list[int]:
    """Group tunnels whose end points are close together.

    Each tunnel is compared, in input order, with the end point of the first
    member (the representative) of every existing cluster, and joins the first
    cluster within ``endpoint_cutoff``; otherwise it starts a new cluster.

    Parameters
    ----------
    tunnels : sequence of Tunnel
    endpoint_cutoff : float, default 5.0
        Maximum end-to-end distance in Å for joining a cluster (inclusive).

    Returns
    -------
    list of int
        Zero-based cluster label for each tunnel, numbered in order of first
        appearance. The result depends on input order.
    """

    labels: list[int] = []
    representatives: list[Tunnel] = []
    for tunnel in tunnels:
        label = None
        for index, representative in enumerate(representatives):
            if distance(tunnel.end, representative.end) <= endpoint_cutoff:
                label = index
                break
        if label is None:
            label = len(representatives)
            representatives.append(tunnel)
        labels.append(label)
    return labels
