"""Immutable data models shared by every CREVICE analysis.

The classes here are small frozen dataclasses. They hold parsed atoms, the
frame they belong to, and the geometric results that the analysis modules
return. Nothing in this module performs an analysis; it defines what the
results mean.

Conventions used throughout
---------------------------
Coordinates
    Cartesian ``(x, y, z)`` tuples in ångström (Å), in the coordinate frame of
    the input file (or of the reference frame after trajectory alignment).
    The alias :data:`Coord` is ``tuple[float, float, float]``.
Clearance
    The signed distance from a point to the nearest atom *surface*, i.e. the
    distance to an atom centre minus that atom's van der Waals radius
    (see :mod:`crevice.radii`). A negative clearance means the point lies
    inside an atom sphere.
Radius
    For channel profiles, the clearance minus the probe radius, clipped at
    zero. For cast/void samples, the clearance itself. Each result class
    below states which definition applies.
Residue labels
    Strings of the form ``"A:ALA42"`` (chain, residue name, residue number,
    optional insertion code); see :attr:`ResidueKey.label`.

All result classes provide ``to_dict()`` for JSON export. The field names in
those dictionaries are the stable output schema; for example the profile
coordinate ``ChannelPoint.t`` is exported as ``"axis_position"``.

These are geometric descriptions of the supplied coordinates. A radius,
volume or contact recorded here is not, on its own, evidence of biological
function, permeation or solvent accessibility.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

#: A Cartesian coordinate or 3-vector ``(x, y, z)``; positions are in Å.
Coord = tuple[float, float, float]


@dataclass(frozen=True, order=True)
class ResidueKey:
    """Hashable identity of one residue.

    Two atoms belong to the same residue exactly when all four fields match.
    Instances are ordered (``order=True``), so sorting a list of keys sorts by
    chain, then residue number, then insertion code, then residue name.

    Attributes
    ----------
    chain_id : str
        Chain identifier as read from the input. May be empty (PDB files without
        chain IDs, or GRO/MD topologies without chain or segment IDs).
    resid : int
        Residue sequence number (author numbering for PDB/mmCIF input).
    icode : str
        PDB insertion code, or ``""`` when absent.
    resname : str
        Residue (component) name, for example ``"ALA"`` or ``"HOH"``.

    Examples
    --------
    >>> from crevice.models import ResidueKey
    >>> ResidueKey("A", 42, "", "ALA").label
    'A:ALA42'
    >>> ResidueKey("", 7, "B", "GLY").label
    '_:GLY7B'
    """
    chain_id: str
    resid: int
    icode: str
    resname: str

    @property
    def label(self) -> str:
        """Human-readable residue label, ``"<chain>:<resname><resid><icode>"``.

        An empty chain ID is written as ``"_"`` and the insertion code is appended
        without a separator. This label is the residue identifier used in every
        CREVICE output (contacts, networks, hydration tables).
        """
        chain = self.chain_id or "_"
        ins = self.icode.strip()
        return f"{chain}:{self.resname}{self.resid}{ins}"


@dataclass(frozen=True)
class Atom:
    """One atom record with coordinates in ångström.

    Attributes
    ----------
    serial : int
        Atom serial number from the input file (PDB ``serial``, mmCIF
        ``_atom_site.id``); for MD frames, the 1-based index within the loaded
        selection.
    name : str
        Atom name, for example ``"CA"`` or ``"OG1"``.
    resname : str
        Residue name.
    chain_id : str
        Chain identifier (may be empty).
    resid : int
        Residue number.
    x, y, z : float
        Cartesian coordinates in Å.
    element : str, default ""
        Chemical element symbol in upper case. The parsers fill this from the
        file or infer it from the atom name (:func:`crevice.radii.infer_element`).
        Hydrogen filtering and van der Waals radii depend on it.
    altloc : str, default ""
        Alternate-location indicator that was retained, or ``""``.
    icode : str, default ""
        Insertion code of the residue.
    occupancy : float or None, default None
        Occupancy from the file, or ``None`` when absent.
    bfactor : float or None, default None
        Isotropic B-factor in Å², or ``None`` when absent.
    hetero : bool, default False
        Whether the atom is treated as a heteroatom. For PDB input this is the
        ``HETATM`` record type; mmCIF input may instead use entity membership
        (see :func:`crevice.mmcif.load_mmcif`); for MD input it means "outside
        the MDAnalysis ``protein`` selection". Analyses use it through
        ``include_hetero`` options.
    """
    serial: int
    name: str
    resname: str
    chain_id: str
    resid: int
    x: float
    y: float
    z: float
    element: str = ""
    altloc: str = ""
    icode: str = ""
    occupancy: float | None = None
    bfactor: float | None = None
    hetero: bool = False

    @property
    def coord(self) -> Coord:
        """Position as an ``(x, y, z)`` tuple in Å."""
        return (self.x, self.y, self.z)

    @property
    def residue_key(self) -> ResidueKey:
        """The :class:`ResidueKey` identifying this atom's residue."""
        return ResidueKey(self.chain_id, self.resid, self.icode, self.resname)


@dataclass(frozen=True)
class Residue:
    """A residue and its atoms, in input order.

    Attributes
    ----------
    key : ResidueKey
        Residue identity.
    atoms : tuple of Atom
        Atoms that share ``key``, in the order they appear in the frame.
    """
    key: ResidueKey
    atoms: tuple[Atom, ...]

    @property
    def label(self) -> str:
        """Residue label, identical to ``key.label``."""
        return self.key.label


@dataclass(frozen=True)
class StructureFrame:
    """One set of atomic coordinates: a static model or one trajectory frame.

    Every analysis in CREVICE takes a ``StructureFrame``. The readers keep
    hydrogens, waters, ions and ligands (they drop only unselected alternate
    locations and, for Gemmi mmCIF input, zero-occupancy atoms; see
    :func:`crevice.parser.load_structure_report`), so hydrogen and heteroatom
    filtering happens later through :meth:`selected_atoms`.

    Attributes
    ----------
    atoms : tuple of Atom
        All atoms in the frame. Must be non-empty.
    source : str, default ""
        Path or description of where the coordinates came from.
    model_index : int, default 0
        Zero-based model number within a multi-model file.
    frame_index : int, default 0
        Zero-based frame number within the source trajectory.

    Raises
    ------
    ValueError
        If ``atoms`` is empty.

    Examples
    --------
    >>> from crevice.models import Atom, StructureFrame
    >>> frame = StructureFrame((
    ...     Atom(1, "N", "GLY", "A", 1, 0.0, 0.0, 0.0, "N"),
    ...     Atom(2, "CA", "GLY", "A", 1, 1.5, 0.0, 0.0, "C"),
    ...     Atom(3, "O", "HOH", "W", 2, 5.0, 0.0, 0.0, "O", hetero=True),
    ... ))
    >>> [r.label for r in frame.residues()]
    ['A:GLY1', 'W:HOH2']
    >>> len(frame.selected_atoms(include_hetero=False))
    2
    >>> frame.centroid()
    (2.1666666666666665, 0.0, 0.0)
    """
    atoms: tuple[Atom, ...]
    source: str = ""
    model_index: int = 0
    frame_index: int = 0

    def __post_init__(self) -> None:
        if not self.atoms:
            raise ValueError("StructureFrame requires at least one atom")

    def residues(self) -> tuple[Residue, ...]:
        """Group the atoms into residues.

        Returns
        -------
        tuple of Residue
            One entry per distinct :class:`ResidueKey`, ordered by the first
            appearance of that residue in :attr:`atoms`. Atoms of a residue that are
            not contiguous in the file are still gathered into one residue.
        """
        grouped: dict[ResidueKey, list[Atom]] = {}
        order: list[ResidueKey] = []
        for atom in self.atoms:
            key = atom.residue_key
            if key not in grouped:
                grouped[key] = []
                order.append(key)
            grouped[key].append(atom)
        return tuple(Residue(key, tuple(grouped[key])) for key in order)

    def selected_atoms(
        self,
        *,
        include_hydrogen: bool = False,
        include_hetero: bool = True,
    ) -> tuple[Atom, ...]:
        """Return the atoms that pass the hydrogen and heteroatom filters.

        Parameters
        ----------
        include_hydrogen : bool, default False
            Keep atoms whose ``element`` is ``"H"`` or deuterium (``"D"``);
            both count as hydrogen.
        include_hetero : bool, default True
            Keep atoms with ``hetero=True`` (ligands, waters, ions, and for MD
            input everything outside the protein selection).

        Returns
        -------
        tuple of Atom
            The retained atoms in input order.

        Raises
        ------
        ValueError
            If no atom survives the filters.
        """
        atoms: list[Atom] = []
        for atom in self.atoms:
            if not include_hetero and atom.hetero:
                continue
            if not include_hydrogen and atom.element.upper() in {"H", "D"}:
                continue
            atoms.append(atom)
        if not atoms:
            raise ValueError("Atom selection is empty")
        return tuple(atoms)

    def bounding_box(self, atoms: Iterable[Atom] | None = None) -> tuple[Coord, Coord]:
        """Axis-aligned bounding box of atom centres.

        Parameters
        ----------
        atoms : iterable of Atom, optional
            Atoms to bound. Defaults to all atoms in the frame.

        Returns
        -------
        tuple of Coord
            ``(min_corner, max_corner)`` in Å. Atom radii are not added.

        Raises
        ------
        ValueError
            If the atom set is empty.
        """
        selected = tuple(atoms) if atoms is not None else self.atoms
        if not selected:
            raise ValueError("Cannot compute a bounding box for an empty atom set")
        xs = [atom.x for atom in selected]
        ys = [atom.y for atom in selected]
        zs = [atom.z for atom in selected]
        return (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))

    def centroid(self, atoms: Iterable[Atom] | None = None) -> Coord:
        """Unweighted mean position of atom centres.

        Parameters
        ----------
        atoms : iterable of Atom, optional
            Atoms to average. Defaults to all atoms in the frame.

        Returns
        -------
        Coord
            Centroid in Å (no mass weighting).

        Raises
        ------
        ValueError
            If the atom set is empty.
        """
        selected = tuple(atoms) if atoms is not None else self.atoms
        if not selected:
            raise ValueError("Cannot compute a centroid for an empty atom set")
        n = float(len(selected))
        return (
            sum(atom.x for atom in selected) / n,
            sum(atom.y for atom in selected) / n,
            sum(atom.z for atom in selected) / n,
        )


@dataclass(frozen=True)
class ChannelPoint:
    """One sample point of a profile, tunnel path or cast.

    The same class is used for several result types; the meaning of
    ``radius`` and ``t`` depends on which analysis produced the point.

    Attributes
    ----------
    index : int
        Zero-based position of the point in its containing sequence.
    position : Coord
        Sample centre ``(x, y, z)`` in Å.
    radius : float
        Radius in Å. For channel profiles
        (:func:`crevice.channels.pore_profile`,
        :func:`crevice.channels.profile_along_path`) this is
        ``max(0, raw_clearance - probe_radius)``: the radius of the largest
        sphere centred here that also leaves room for the probe. For void and
        rolling-probe cast samples and tunnel grid nodes it equals
        ``raw_clearance``.
    raw_clearance : float
        Signed distance in Å from ``position`` to the nearest atom van der
        Waals surface (centre distance minus atomic radius). Negative inside an
        atom sphere.
    t : float
        Coordinate along the result's reference direction, in Å unless stated.
        For axial profiles it is the signed axial position relative to
        ``PoreProfile.axis_origin``; for explicit paths it is the cumulative
        path length from the first point; for
        :func:`crevice.tunnels.find_tunnels` paths it is also the cumulative
        path length from the start (the step number is :attr:`index`); for
        cast samples it is the projection on the cast's axis or grid
        direction. Exported as ``"axis_position"``.
    nearest_atom_serial : int or None, default None
        Serial of the atom whose surface is nearest to ``position``.
    nearest_residue : str or None, default None
        Label (:attr:`ResidueKey.label`) of that atom's residue.
    """
    index: int
    position: Coord
    radius: float
    raw_clearance: float
    t: float
    nearest_atom_serial: int | None = None
    nearest_residue: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary.

        Returns
        -------
        dict
            Keys ``index``, ``x``, ``y``, ``z``, ``radius``, ``raw_clearance``,
            ``axis_position`` (the value of :attr:`t`), ``nearest_atom_serial`` and
            ``nearest_residue``.
        """
        return {
            "index": self.index,
            "x": self.position[0],
            "y": self.position[1],
            "z": self.position[2],
            "radius": self.radius,
            "raw_clearance": self.raw_clearance,
            "axis_position": self.t,
            "nearest_atom_serial": self.nearest_atom_serial,
            "nearest_residue": self.nearest_residue,
        }


@dataclass(frozen=True)
class PoreProfile:
    """A radius profile sampled along a channel axis or an explicit path.

    Attributes
    ----------
    axis_origin : Coord
        Reference point of the axis in Å. The ``t`` of each point is measured
        from here.
    axis_direction : Coord
        Unit direction of the axis. For explicit paths and tunnels it is the
        unit vector from the first to the last point.
    points : tuple of ChannelPoint
        Samples in profile order. ``radius`` is clearance minus
        :attr:`probe_radius`, clipped at zero.
    probe_radius : float, default 0.0
        Probe radius in Å that was subtracted from the raw clearance.
    method : str, default "axial-refined"
        How the profile was obtained: ``"axial-connected"`` (connected channel
        sections, the default of :func:`crevice.channels.pore_profile`),
        ``"axial-refined"`` (fixed-axis scan, ``search_radius=0``),
        ``"explicit-path"`` (:func:`crevice.channels.profile_along_path`) or
        ``"tunnel-path"`` (:func:`crevice.tunnels.tunnel_profile`).
    metadata : dict
        Settings and diagnostics recorded by the producing function (atom
        selection, status, sampling). Contents vary by method.

    Notes
    -----
    Summary properties (:attr:`min_radius`, :attr:`mean_radius`, ...) are
    computed from the stored samples only. They are sampling estimates: a
    narrower constriction between two samples is not detected.
    """
    axis_origin: Coord
    axis_direction: Coord
    points: tuple[ChannelPoint, ...]
    probe_radius: float = 0.0
    method: str = "axial-refined"
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def bottleneck(self) -> ChannelPoint:
        """The sample with the smallest ``radius``.

        Ties are resolved in favour of the earliest sample. Raises ``ValueError``
        for a profile without points.
        """
        return min(self.points, key=lambda point: point.radius)

    @property
    def min_radius(self) -> float:
        """Smallest sampled radius in Å (the bottleneck radius)."""
        return self.bottleneck.radius

    @property
    def max_radius(self) -> float:
        """Largest sampled radius in Å."""
        return max(point.radius for point in self.points)

    @property
    def mean_radius(self) -> float:
        """Arithmetic mean of the sampled radii in Å (not length weighted)."""
        return sum(point.radius for point in self.points) / len(self.points)

    @property
    def length(self) -> float:
        """Polyline length through the sample positions, in Å.

        Returns 0.0 when there are fewer than two points. This follows the refined
        centres, so it can exceed the axial extent of the profile.
        """
        from .geometry import distance

        if len(self.points) < 2:
            return 0.0
        return sum(
            distance(left.position, right.position)
            for left, right in zip(self.points, self.points[1:])
        )

    @property
    def volume_estimate(self) -> float:
        """Volume enclosed by circular cross-sections along the sampled path, in Å³.

        Each segment between consecutive samples contributes
        ``pi * (r_i**2 + r_j**2) / 2 * |p_j - p_i|`` (trapezoidal integration of the
        cross-sectional area ``pi * r**2``). Returns 0.0 for fewer than two points.

        Notes
        -----
        This assumes circular cross-sections of the probe-corrected radius, so it is
        a coarse estimate for non-circular channels. For a measured volume use a
        cast (:func:`crevice.rolling.rolling_probe_cast`).
        """
        if len(self.points) < 2:
            return 0.0
        import math

        from .geometry import distance

        # Trapezoidal integration of sampled cross-sectional areas.
        return sum(math.pi * (left.radius**2 + right.radius**2) * 0.5
                   * distance(left.position, right.position)
                   for left, right in zip(self.points, self.points[1:]))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary including summaries and all points.

        Returns
        -------
        dict
            ``method``, axis, ``probe_radius``, ``point_count``, ``min_radius``,
            ``max_radius``, ``mean_radius``, ``length``, ``volume_estimate``,
            ``bottleneck`` (a point dictionary), ``metadata`` and ``points``.
        """
        bottleneck = self.bottleneck
        return {
            "method": self.method,
            "axis_origin": self.axis_origin,
            "axis_direction": self.axis_direction,
            "probe_radius": self.probe_radius,
            "point_count": len(self.points),
            "min_radius": self.min_radius,
            "max_radius": self.max_radius,
            "mean_radius": self.mean_radius,
            "length": self.length,
            "volume_estimate": self.volume_estimate,
            "bottleneck": bottleneck.to_dict(),
            "metadata": self.metadata,
            "points": [point.to_dict() for point in self.points],
        }


@dataclass(frozen=True)
class Cavity:
    """Summary of one cavity found by :func:`crevice.cavities.detect_cavities`.

    Attributes
    ----------
    id : int
        One-based identifier in the order returned (largest volume first).
    center : Coord
        Clearance-weighted centroid of the cavity's grid points, in Å.
    radius : float
        Largest atom-surface clearance among the cavity's grid points, in Å:
        the radius of the largest atom-free sphere centred on a grid point.
    volume : float
        ``grid_points * spacing**3`` in Å³: the volume of the grid cells whose
        centres satisfy the clearance threshold. Depends on grid spacing.
    grid_points : int
        Number of grid points in the cavity.
    nearest_residues : tuple of str, default ()
        Up to ten residue labels within the shell radius of ``center``, closest
        first.
    score : float, default 0.0
        Ranking value copied from the cast component. For
        :func:`crevice.cavities.detect_cavities` it equals ``volume``. Not a
        probability or energy.
    """
    id: int
    center: Coord
    radius: float
    volume: float
    grid_points: int
    nearest_residues: tuple[str, ...] = ()
    score: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary of all fields."""
        return {
            "id": self.id,
            "center": self.center,
            "radius": self.radius,
            "volume": self.volume,
            "grid_points": self.grid_points,
            "nearest_residues": list(self.nearest_residues),
            "score": self.score,
        }


@dataclass(frozen=True)
class Tunnel:
    """A grid path from a start point to the edge of the search box.

    Produced by :func:`crevice.tunnels.find_tunnels`.

    Attributes
    ----------
    id : int
        One-based identifier in the order found (widest bottleneck first).
    start : Coord
        First point of the path in Å (the requested start, or the molecule
        centroid by default).
    end : Coord
        Grid node on the box boundary where the path ends, in Å.
    points : tuple of ChannelPoint
        Path samples. ``radius`` and ``raw_clearance`` are the atom-surface
        clearance at each node; ``t`` is the cumulative path length in Å from
        the start and ``index`` the step number.
    length : float
        Polyline length of the path in Å.
    bottleneck_radius : float
        Minimum clearance in Å over the path, taking the whole connecting
        segments into account (not only the nodes).
    mean_radius : float
        Arithmetic mean of node clearances in Å.
    score : float
        ``bottleneck_radius / (1 + 0.01 * length)``: a ranking heuristic that
        favours wide, short paths.
    nearest_residues : tuple of str, default ()
        Up to ten residue labels lining the path.
    metadata : dict
        Grid spacing, seed attachment, per-segment clearances and estimator
        definitions.
    """
    id: int
    start: Coord
    end: Coord
    points: tuple[ChannelPoint, ...]
    length: float
    bottleneck_radius: float
    mean_radius: float
    score: float
    nearest_residues: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary of all fields, including points."""
        return {
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "bottleneck_radius": self.bottleneck_radius,
            "mean_radius": self.mean_radius,
            "score": self.score,
            "nearest_residues": list(self.nearest_residues),
            "points": [point.to_dict() for point in self.points],
            "metadata": self.metadata,
        }


@dataclass(frozen=True)
class VoidComponent:
    """One connected region of a cast (void, cavity or channel fill).

    Attributes
    ----------
    id : int
        One-based identifier within its :class:`VoidCast`.
    kind : str
        Topological or selection class, for example ``"buried_cavity"``,
        ``"surface_connected"``, ``"through_channel"``,
        ``"centerline_channel"``, ``"rolling_probe_void"`` or
        ``"dominant_cavity_region"``. Values are geometric labels only.
    center : Coord
        Representative centre in Å (a clearance-weighted centroid for grid
        casts, the plain centroid for rolling-probe casts).
    volume : float
        Sampled volume in Å³: number of samples times ``spacing**3``.
    points : tuple of ChannelPoint
        Every measurement sample of the component. ``radius`` equals the
        sample's atom-surface clearance.
    boundary_faces : tuple of str, default ()
        Faces of the grid box touched by the component (``"x-"``, ``"x+"``,
        ...). Empty for rolling-probe casts.
    nearest_residues : tuple of str, default ()
        Residue labels near the component.
    score : float, default 0.0
        Ranking heuristic of the producing method (not a probability).
    metadata : dict
        Method-specific diagnostics.
    """
    id: int
    kind: str
    center: Coord
    volume: float
    points: tuple[ChannelPoint, ...]
    boundary_faces: tuple[str, ...] = ()
    nearest_residues: tuple[str, ...] = ()
    score: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def point_count(self) -> int:
        """Number of samples in the component."""
        return len(self.points)

    @property
    def min_radius(self) -> float:
        """Smallest sample radius in Å, or 0.0 when there are no points."""
        return min((point.radius for point in self.points), default=0.0)

    @property
    def mean_radius(self) -> float:
        """Mean sample radius in Å, or 0.0 when there are no points."""
        if not self.points:
            return 0.0
        return sum(point.radius for point in self.points) / len(self.points)

    @property
    def max_radius(self) -> float:
        """Largest sample radius in Å, or 0.0 when there are no points."""
        return max((point.radius for point in self.points), default=0.0)

    def to_dict(self, *, include_points: bool = False) -> dict[str, Any]:
        """Return a JSON-friendly dictionary.

        Parameters
        ----------
        include_points : bool, default False
            Also include every sample under ``"points"``. Off by default because
            casts can contain many thousands of samples.

        Returns
        -------
        dict
            Component fields plus ``point_count`` and the radius summaries.
        """
        data: dict[str, Any] = {
            "id": self.id,
            "kind": self.kind,
            "center": self.center,
            "volume": self.volume,
            "point_count": self.point_count,
            "min_radius": self.min_radius,
            "mean_radius": self.mean_radius,
            "max_radius": self.max_radius,
            "boundary_faces": list(self.boundary_faces),
            "nearest_residues": list(self.nearest_residues),
            "score": self.score,
            "metadata": self.metadata,
        }
        if include_points:
            data["points"] = [point.to_dict() for point in self.points]
        return data


@dataclass(frozen=True)
class VoidCast:
    """Result of a cast: selected components plus display samples.

    Produced by :func:`crevice.voids.void_cast`,
    :func:`crevice.rolling.rolling_probe_cast` and
    :func:`crevice.regional_volume.regional_probe_cast`.

    Attributes
    ----------
    components : tuple of VoidComponent
        Retained components. Their ``points`` are the full measurement samples.
    points : tuple of ChannelPoint
        Samples for display and export. They may be an evenly thinned subset of
        the component samples (see ``metadata["display_subsampled"]``);
        thinning never changes the component volumes.
    spacing : float
        Grid spacing in Å actually used (after any automatic coarsening).
    min_radius : float
        Clearance threshold in Å that every grid-cast sample satisfies. Rolling
        and regional casts record 0.0 here and store their probe radius in
        ``metadata["probe_radius_A"]``.
    mode : str
        Cast mode (``"auto"``, ``"channel"``, ``"cavity"``, ``"all"`` or
        ``"rolling"``).
    metadata : dict
        Parameters, grid definition and selection decisions of the producing
        method.
    """
    components: tuple[VoidComponent, ...]
    points: tuple[ChannelPoint, ...]
    spacing: float
    min_radius: float
    mode: str
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def point_count(self) -> int:
        """Number of display samples in :attr:`points`."""
        return len(self.points)

    @property
    def component_count(self) -> int:
        """Number of retained components."""
        return len(self.components)

    @property
    def total_volume(self) -> float:
        """Sum of component volumes in Å³ (independent of display thinning)."""
        return sum(component.volume for component in self.components)

    @property
    def bottleneck(self) -> ChannelPoint | None:
        """Display sample with the smallest radius, or ``None`` when there are no samples."""
        if not self.points:
            return None
        return min(self.points, key=lambda point: point.radius)

    @property
    def min_clearance(self) -> float:
        """Smallest radius among the display samples in Å (0.0 when empty)."""
        return min((point.radius for point in self.points), default=0.0)

    @property
    def mean_clearance(self) -> float:
        """Mean radius of the display samples in Å (0.0 when empty)."""
        if not self.points:
            return 0.0
        return sum(point.radius for point in self.points) / len(self.points)

    @property
    def max_clearance(self) -> float:
        """Largest radius among the display samples in Å (0.0 when empty)."""
        return max((point.radius for point in self.points), default=0.0)

    def to_dict(self, *, include_points: bool = False) -> dict[str, Any]:
        """Return a JSON-friendly dictionary.

        The ``"method"`` key is always ``"clearance-grid-fill"``; the producing
        algorithm is recorded in ``metadata``.

        Parameters
        ----------
        include_points : bool, default False
            Also include the display samples under ``"points"``.

        Returns
        -------
        dict
            Cast settings, summaries, ``metadata`` and component dictionaries
            (without their points).
        """
        data: dict[str, Any] = {
            "method": "clearance-grid-fill",
            "mode": self.mode,
            "spacing": self.spacing,
            "min_radius": self.min_radius,
            "component_count": self.component_count,
            "point_count": self.point_count,
            "total_volume": self.total_volume,
            "min_clearance": self.min_clearance,
            "mean_clearance": self.mean_clearance,
            "max_clearance": self.max_clearance,
            "bottleneck": self.bottleneck.to_dict() if self.bottleneck else None,
            "metadata": self.metadata,
            "components": [component.to_dict() for component in self.components],
        }
        if include_points:
            data["points"] = [point.to_dict() for point in self.points]
        return data


@dataclass(frozen=True)
class ResidueContact:
    """Geometric contact between one residue and a sampled region.

    Produced by :func:`crevice.analysis.annotate_residues`.

    Attributes
    ----------
    residue : str
        Residue label.
    min_distance : float
        Smallest gap in Å between the residue's atom surfaces and the region
        surface (for profile points, the sphere of radius ``raw_clearance``
        around each sample; for bare coordinates, the point itself), clipped at
        zero.
    mean_distance : float
        Mean of those gaps over all atom/point pairs within the cutoff, in Å.
    atom_count : int
        Number of the residue's atoms in the analysed selection.
    point_count : int
        Number of region samples within the cutoff of the residue.
    role : str
        ``"bottleneck"`` (``min_distance <= 1 Å`` and within two samples of the
        bottleneck index), ``"lining"`` (``min_distance <= 1 Å``),
        ``"bottleneck-nearby"`` (near the bottleneck but further than 1 Å) or
        ``"nearby"``.
    properties : tuple of str, default ()
        Residue classes from :func:`crevice.radii.residue_properties`.
    influence_score : float, default 0.0
        Ranking heuristic ``1/(0.25 + min_distance) + 0.15*log(1 + point_count)``
        plus 1.0 for ``"bottleneck"`` and ``"bottleneck-nearby"`` and 0.25 for ``"lining"``. Describes
        geometric proximity only.
    """
    residue: str
    min_distance: float
    mean_distance: float
    atom_count: int
    point_count: int
    role: str
    properties: tuple[str, ...] = ()
    influence_score: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary of all fields."""
        return {
            "residue": self.residue,
            "min_distance": self.min_distance,
            "mean_distance": self.mean_distance,
            "atom_count": self.atom_count,
            "point_count": self.point_count,
            "role": self.role,
            "properties": list(self.properties),
            "influence_score": self.influence_score,
        }

@dataclass(frozen=True)
class NetworkNode:
    """Node of a residue-interaction network.

    Attributes
    ----------
    id : str
        Unique node identifier (the residue label, or ``"region:..."`` for a
        channel/cavity pseudo-node).
    kind : str
        ``"residue"`` or ``"region"``.
    label : str
        Display label.
    residue : str or None, default None
        Residue label for residue nodes.
    properties : tuple of str, default ()
        Residue classes (see :func:`crevice.radii.residue_properties`).
    metadata : dict
        For residue nodes: ``resname``, ``chain_id``, ``resid``, ``icode``,
        user ``group`` and ``position`` (centroid of the selected atoms, Å).
    """
    id: str
    kind: str
    label: str
    residue: str | None = None
    properties: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary of all fields."""
        return {
            "id": self.id,
            "kind": self.kind,
            "label": self.label,
            "residue": self.residue,
            "properties": list(self.properties),
            "metadata": self.metadata,
        }


@dataclass(frozen=True)
class NetworkEdge:
    """Undirected edge of a residue-interaction network.

    Attributes
    ----------
    source, target : str
        Node identifiers. For residue-residue edges ``source < target``.
    interaction : str
        Geometric/chemical class, for example ``"salt_bridge_candidate"``,
        ``"hydrophobic_residue_contact"`` or, for region edges,
        ``"region_lining"``. See :func:`crevice.networks.build_residue_network`.
    distance : float
        Distance in Å under the network's metric (surface gap or atom-centre
        distance; see ``ResidueInteractionNetwork.metadata["distance_metric"]``).
    weight : float, default 1.0
        Edge weight; ``1/(1 + distance)`` for residue edges (surface gaps below
        zero count as zero), the contact influence score for region edges.
    metadata : dict
        Closest atom pair and other evidence.
    """
    source: str
    target: str
    interaction: str
    distance: float
    weight: float = 1.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary of all fields."""
        return {
            "source": self.source,
            "target": self.target,
            "interaction": self.interaction,
            "distance": self.distance,
            "weight": self.weight,
            "metadata": self.metadata,
        }


@dataclass(frozen=True)
class ResidueInteractionNetwork:
    """A set of residue (and optional region) nodes joined by contact edges.

    Attributes
    ----------
    nodes : tuple of NetworkNode
    edges : tuple of NetworkEdge
    metadata : dict
        Construction settings: ``cutoff`` (Å), ``distance_metric``,
        atom-selection flags and rule identifiers.
    """
    nodes: tuple[NetworkNode, ...]
    edges: tuple[NetworkEdge, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def degree(self) -> dict[str, int]:
        """Number of edges incident on each node.

        Returns
        -------
        dict of str to int
            Every node appears, with 0 for isolated nodes. An edge whose endpoint is
            not in :attr:`nodes` still increments that endpoint's count.
        """
        counts = {node.id: 0 for node in self.nodes}
        for edge in self.edges:
            counts[edge.source] = counts.get(edge.source, 0) + 1
            counts[edge.target] = counts.get(edge.target, 0) + 1
        return counts

    def neighbors(self, node_id: str) -> tuple[str, ...]:
        """Identifiers of nodes that share an edge with ``node_id``.

        Parameters
        ----------
        node_id : str
            Node identifier.

        Returns
        -------
        tuple of str
            Sorted, without duplicates. Empty if the node has no edges or does not
            exist.
        """
        neighbors: set[str] = set()
        for edge in self.edges:
            if edge.source == node_id:
                neighbors.add(edge.target)
            elif edge.target == node_id:
                neighbors.add(edge.source)
        return tuple(sorted(neighbors))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary with counts, metadata, nodes and edges."""
        return {
            "node_count": len(self.nodes),
            "edge_count": len(self.edges),
            "metadata": self.metadata,
            "nodes": [node.to_dict() for node in self.nodes],
            "edges": [edge.to_dict() for edge in self.edges],
        }


@dataclass(frozen=True)
class FrameAnalysis:
    """All results computed for one trajectory frame.

    Produced by :func:`crevice.trajectory.analyze_trajectory`.

    Attributes
    ----------
    frame_index : int
        Frame number in the source trajectory.
    time : float or None, default None
        Frame time as reported by the reader (ps for MDAnalysis input), or
        ``None`` when unknown.
    profile : PoreProfile or None, default None
        Channel profile, or ``None`` when not requested or unresolved.
    cavities : tuple of Cavity, default ()
    tunnels : tuple of Tunnel, default ()
    residue_contacts : tuple of ResidueContact, default ()
    network : ResidueInteractionNetwork or None, default None
    features : dict of str to float
        Numeric features (see :func:`crevice.ml.profile_features`).
    metadata : dict
        Alignment report, profile status and reason, constriction evidence and
        network counts.
    """
    frame_index: int
    time: float | None = None
    profile: PoreProfile | None = None
    cavities: tuple[Cavity, ...] = ()
    tunnels: tuple[Tunnel, ...] = ()
    residue_contacts: tuple[ResidueContact, ...] = ()
    network: ResidueInteractionNetwork | None = None
    features: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary with nested result dictionaries."""
        return {
            "frame_index": self.frame_index,
            "time": self.time,
            "profile": self.profile.to_dict() if self.profile else None,
            "cavities": [cavity.to_dict() for cavity in self.cavities],
            "tunnels": [tunnel.to_dict() for tunnel in self.tunnels],
            "residue_contacts": [contact.to_dict() for contact in self.residue_contacts],
            "network": self.network.to_dict() if self.network else None,
            "features": self.features,
            "metadata": self.metadata,
        }


@dataclass(frozen=True)
class TrajectoryAnalysis:
    """Per-frame results of a multi-frame analysis.

    Attributes
    ----------
    frames : tuple of FrameAnalysis
        One entry per analysed frame, in order.
    analyses : tuple of str
        Names of the analyses that were run.
    metadata : dict
        Trajectory-level settings and summaries (stride, alignment, reference
        profile status, constriction residue frequencies).
    """
    frames: tuple[FrameAnalysis, ...]
    analyses: tuple[str, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly dictionary including every frame."""
        return {
            "frame_count": len(self.frames),
            "analyses": list(self.analyses),
            "metadata": self.metadata,
            "frames": [frame.to_dict() for frame in self.frames],
        }
