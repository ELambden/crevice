"""Shared helpers for the non-geometric validation scripts in ``benchmarks/nongeom``.

Three kinds of helper live here:

* builders of synthetic MD systems (an MDAnalysis ``Universe`` in memory with
  a fixed protein ring for alignment and one residue per scripted molecule);
* simple stochastic models with known answers (independent and single-file
  random walkers on a periodic line, drift-diffusion of ions);
* an **independent reference implementation** of permeation counting and of
  the collective coordinate ``n(t)``, written from the published definitions
  (Zhu, Tajkhorshid and Schulten, Biophys. J. 86, 50-57, 2004; Phys. Rev.
  Lett. 93, 224501, 2004) using only NumPy and MDAnalysis coordinates. It does
  not import or call anything from ``crevice``; it exists so that CREVICE's
  counts can be compared with a second program written separately.

Nothing here is a biological validation: the systems are constructed so that
the right answer is known by construction or analytically.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

NONGEOM = Path(__file__).resolve().parent
BENCH = NONGEOM.parent
REPO = BENCH.parent
RESULTS = BENCH / "results" / "nongeom"
sys.path.insert(0, str(BENCH))
from common import machine_record, write_csv, write_json  # noqa: E402,F401

#: Siegmund's constant: a Gaussian random walk monitored at discrete steps sees
#: an absorbing barrier shifted outward by beta * sigma (sigma = step SD);
#: beta = -zeta(1/2) / sqrt(2 pi) (Siegmund, Adv. Appl. Probab. 11, 701, 1979).
SIEGMUND_BETA = 0.5825971579390106

STATE_BELOW, STATE_INSIDE, STATE_ABOVE, STATE_OUT = -1, 0, 1, 2


def ensure_src_on_path() -> None:
    """Make ``src/crevice`` importable for the scripts (the reference code never needs it)."""
    src = str(REPO / "src")
    if src not in sys.path:
        sys.path.insert(0, src)


# ---------------------------------------------------------------------------
# Synthetic systems
# ---------------------------------------------------------------------------

def ring_atoms(radius=6.0, count=6):
    """Positions of a fixed 'protein' ring of CA atoms about the origin (used for alignment)."""
    return [(radius * math.cos(a), radius * math.sin(a), 0.0)
            for a in np.linspace(0, 2 * math.pi, count, endpoint=False)]


def memory_universe(molecule_resnames, positions, *, box, dt=1.0, protein=None, atom_names=None,
                    centre=None, wrap=True):
    """MDAnalysis Universe with a fixed CA ring plus one single-atom residue per molecule.

    Parameters
    ----------
    molecule_resnames : sequence of str
        Residue name of each molecule (``TIP3`` water, ``SOD``/``CLA`` ions, ...).
    positions : array, shape (frames, molecules, 3)
        Positions relative to ``centre`` (default the box centre), Å.
    box : sequence of 3 float
        Orthorhombic box lengths, Å.
    dt : float
        Frame interval, ps.
    protein : sequence of (x, y, z), optional
        Protein CA positions relative to ``centre`` (default :func:`ring_atoms`).
    atom_names : sequence of str, optional
        Atom name per molecule (default ``OH2`` for water, else the residue name).
    wrap : bool
        Wrap all coordinates into the primary box (as an MD engine writes them).
    """
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    positions = np.asarray(positions, dtype=float)
    frames, n_mol, _ = positions.shape
    protein = ring_atoms() if protein is None else protein
    n_prot = len(protein)
    names = ["CA"] * n_prot
    if atom_names is None:
        atom_names = ["OH2" if r in ("TIP3", "SOL", "WAT") else r for r in molecule_resnames]
    names += list(atom_names)
    resnames = ["ALA"] * n_prot + list(molecule_resnames)
    n_atoms = n_prot + n_mol
    u = mda.Universe.empty(n_atoms, n_residues=n_atoms, atom_resindex=np.arange(n_atoms), trajectory=False)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", np.arange(1, n_atoms + 1))
    u.add_TopologyAttr("segids", ["SYS"])
    box = np.asarray(box, dtype=float)
    c = box / 2 if centre is None else np.asarray(centre, dtype=float)
    xyz = np.empty((frames, n_atoms, 3), dtype=np.float32)
    xyz[:, :n_prot] = (np.asarray(protein, dtype=float) + c)[None]
    xyz[:, n_prot:] = positions + c
    if wrap:
        xyz = np.mod(xyz, box.astype(np.float32))
    u.load_new(xyz, format=MemoryReader, dimensions=np.array([*box, 90, 90, 90], np.float32), dt=dt)
    return u


def write_md_files(universe, directory, stem="system"):
    """Write a universe as PDB topology + XTC trajectory (for CLI-level runs)."""
    import MDAnalysis as mda
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    top, traj = directory / f"{stem}.pdb", directory / f"{stem}.xtc"
    universe.trajectory[0]
    universe.atoms.write(str(top))
    with mda.Writer(str(traj), universe.atoms.n_atoms) as writer:
        for _ in universe.trajectory:
            writer.write(universe.atoms)
    return top, traj


# ---------------------------------------------------------------------------
# Stochastic models with known answers
# ---------------------------------------------------------------------------

def ring_walkers(n_particles, frames, *, D, dt, box_z, drift=0.0, rng, start=None):
    """Unwrapped 1D Brownian positions on a periodic line of length ``box_z``.

    ``x(t+dt) = x(t) + drift*dt + N(0, 2 D dt)``; returns shape (frames, n).
    Positions start uniformly over the period unless ``start`` is given.
    """
    x0 = rng.uniform(-box_z / 2, box_z / 2, n_particles) if start is None else np.asarray(start, float)
    steps = rng.normal(drift * dt, math.sqrt(2 * D * dt), size=(frames - 1, n_particles))
    return np.vstack([x0[None], x0[None] + np.cumsum(steps, axis=0)])


def single_file_relabel(x, box_z):
    """Point particles with hard-core reflection on a ring, from independent walkers.

    For point Brownian particles, reflecting collisions are equivalent to
    exchanging labels of independent walkers. At every frame the labelled
    positions are the cyclically ordered values of ``x mod box_z``, unwrapped
    with the offset that conserves the total displacement (sum of unwrapped
    positions), so the collective coordinate is unchanged by the relabelling
    while tagged particles can no longer pass each other.

    Parameters
    ----------
    x : array, shape (frames, n)
        Unwrapped independent-walker positions.
    box_z : float
        Period.

    Returns
    -------
    numpy.ndarray, shape (frames, n)
        Unwrapped single-file positions, ordered y_0 < y_1 < ... < y_0 + box_z.
    """
    x = np.asarray(x, dtype=float)
    frames, n = x.shape
    out = np.empty_like(x)
    total = x.sum(axis=1)
    for t in range(frames):
        base = np.sort(np.mod(x[t], box_z))
        # s_j = base[j mod n] + floor(j / n) * box_z; choose m with sum_{k} s_{k+m} = total.
        m = int(round((total[t] - base.sum()) / box_z))
        j = np.arange(n) + m
        out[t] = base[np.mod(j, n)] + np.floor_divide(j, n) * box_z
    # The sum-conserving offset keeps labels continuous: when a particle crosses the
    # period boundary, its base index and m change together, so no label jumps.
    return out


def ring_n_msd(tau, n_particles, D, L, period, kmax=4000):
    """Exact MSD of Zhu's n(t) for independent (or single-file point) walkers on a ring.

    Each walker contributes g(u) = u / period + f(u), with f periodic, so
    MSD(tau) = N [2 D tau / period^2 + 2 sum_{k != 0} |f_k|^2 (1 - exp(-D q_k^2 tau))],
    q_k = 2 pi k / period, |f_k|^2 = 4 sin^2(q_k L / 2) / (q_k^4 period^2 L^2). The
    cross term vanishes for uniformly distributed walkers. Short-time slope:
    2 N D / (L period); long-time slope: 2 N D / period^2. Hard-core relabelling
    does not change n, so the result also holds for point particles in single file.
    """
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    k = np.arange(1, kmax + 1)
    q = 2 * np.pi * k / period
    f2 = 4 * np.sin(q * L / 2) ** 2 / (q ** 4 * period ** 2 * L ** 2)
    bounded = 2 * 2 * np.sum(f2[None, :] * (1 - np.exp(-D * q[None, :] ** 2 * tau[:, None])), axis=1)
    return n_particles * (2 * D * tau / period ** 2 + bounded)


def expected_window_fit(lag_lo, lag_hi, dt, n_particles, D, L, period):
    """D_n that a least-squares fit a + 2 D tau over lags lag_lo..lag_hi recovers from the exact MSD."""
    lags = np.arange(lag_lo, lag_hi + 1) * dt
    slope = np.polyfit(lags, ring_n_msd(lags, n_particles, D, L, period), 1)[0]
    return slope / 2


def ou_noise(frames, *, tau_frames, sd, rng):
    """Stationary Ornstein-Uhlenbeck (AR(1)) series with correlation time ``tau_frames`` and SD ``sd``."""
    if tau_frames <= 0:
        return rng.normal(0, sd, frames)
    phi = math.exp(-1.0 / tau_frames)
    e = rng.normal(0, sd * math.sqrt(1 - phi ** 2), frames)
    y = np.empty(frames)
    y[0] = rng.normal(0, sd)
    for t in range(1, frames):
        y[t] = phi * y[t - 1] + e[t]
    return y


def ar1(n, phi, rng, sigma=1.0):
    """Stationary AR(1) series ``x_t = phi x_(t-1) + e_t`` with marginal SD ``sigma``."""
    e = rng.normal(0, sigma * math.sqrt(1 - phi ** 2), n)
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x


def ar1_true_tau(phi):
    """Integrated autocorrelation time of AR(1): (1 + phi) / (1 - phi)."""
    return (1 + phi) / (1 - phi)


def ar1_true_se(n, phi, sigma=1.0):
    """Exact SD of the mean of n stationary AR(1) samples."""
    k = np.arange(1, n)
    return sigma * math.sqrt((1 + 2 * np.sum((1 - k / n) * phi ** k)) / n)


# ---------------------------------------------------------------------------
# Independent reference permeation counter (no crevice imports)
# ---------------------------------------------------------------------------

def reference_axial_lateral(universe, atom_indices, centre, *, axis_index=2, stride=1, start=0, stop=None):
    """Unwrapped axial coordinate and lateral distance of molecules over a trajectory.

    Orthorhombic minimum image about ``centre`` for the first frame, then the
    frame-to-frame axial displacement is minimum-imaged and accumulated
    (unwrapped). No alignment (the synthetic systems do not move as a whole).

    Returns
    -------
    s : array (frames, n)   unwrapped axial coordinate (absolute, Å)
    s_img : array (frames, n) axial coordinate of the image nearest ``centre``
    lateral : array (frames, n) lateral distance from the axis line through ``centre``
    times : array (frames,)
    """
    centre = np.asarray(centre, dtype=float)
    lat_axes = [k for k in range(3) if k != axis_index]
    atoms = universe.atoms[np.asarray(atom_indices)]
    s_list, img_list, lat_list, times = [], [], [], []
    prev_img = None
    for ts in universe.trajectory[start:stop:stride]:
        box = np.asarray(ts.dimensions[:3], dtype=float)
        d = atoms.positions.astype(float) - centre
        d -= box * np.round(d / box)
        img = centre[axis_index] + d[:, axis_index]
        if prev_img is None:
            s = img.copy()
        else:
            step = img - prev_img
            L = box[axis_index]
            step -= L * np.round(step / L)
            s = s_list[-1] + step
        prev_img = img
        s_list.append(s)
        img_list.append(img)
        lat_list.append(np.hypot(d[:, lat_axes[0]], d[:, lat_axes[1]]))
        times.append(float(ts.time))
    return np.asarray(s_list), np.asarray(img_list), np.asarray(lat_list), np.asarray(times)


def reference_states(s_img, lateral, zmin, zmax, radius):
    """Per-frame state below/inside/above/out from the image coordinate and lateral distance."""
    state = np.where(s_img < zmin, STATE_BELOW, np.where(s_img > zmax, STATE_ABOVE, STATE_OUT))
    slab = state == STATE_OUT
    ok = np.ones_like(slab) if radius is None else lateral <= radius
    state[slab & ok] = STATE_INSIDE
    return state


def _runs(seq):
    """Run-length encode a 1D sequence: list of (value, first index, last index)."""
    out = []
    start = 0
    for i in range(1, len(seq) + 1):
        if i == len(seq) or seq[i] != seq[start]:
            out.append((int(seq[start]), start, i - 1))
            start = i
    return out


def reference_events(states, s_unwrapped, zmin, zmax, *, convention="strict", s_image=None):
    """Permeation events from per-frame states, written from the definition.

    Conventions
    -----------
    ``strict``
        Zhu et al. (2004): bulk side -> inside the lumen without interruption ->
        opposite bulk side. A frame outside the lumen laterally ("out") ends the
        visit and forgets the side.
    ``compartment`` (``side_memory`` is an alias)
        The molecule remembers the bulk side it last occupied. A passage is a
        run of frames within the axial extent between two bulk frames; when
        it ends on the side opposite to its origin it is an event, attributed
        to the channel ("complete") if at least one frame of the passage was
        inside the lumen and listed as "outside_lumen" otherwise. This is the
        compartment rule of GROMACS computational electrophysiology.

    A whole-extent jump between two consecutive frames is an *unresolved*
    crossing (counted separately, never as an event): the previous image
    coordinate (``s_image``, default ``s_unwrapped``) plus the unwrapped step
    lies beyond the far end of the extent. Returns a list of dicts with molecule, direction, entry and
    exit frame, status.
    """
    frames, n = states.shape
    events = []
    if convention == "side_memory":
        convention = "compartment"
    if convention not in ("strict", "compartment"):
        raise ValueError(convention)
    for i in range(n):
        seq = states[:, i].copy()
        if convention == "compartment":
            bulk = (seq == STATE_BELOW) | (seq == STATE_ABOVE)
            last_bulk = None          # (frame, side) of the latest bulk frame
            for f in range(frames):
                if not bulk[f]:
                    continue
                if last_bulk is not None and f > last_bulk[0] + 1 and seq[f] == -last_bulk[1]:
                    passage = seq[last_bulk[0] + 1:f]
                    events.append({"molecule": i, "direction": "up" if last_bulk[1] == STATE_BELOW else "down",
                                   "entry_frame": last_bulk[0] + 1, "exit_frame": f,
                                   "status": "complete" if np.any(passage == STATE_INSIDE) else "outside_lumen"})
                last_bulk = (f, int(seq[f]))
        else:
            runs = _runs(seq) if len(seq) else []
            for r in range(1, len(runs) - 1):
                value, a, b = runs[r]
                if value != STATE_INSIDE:
                    continue
                before, after = runs[r - 1][0], runs[r + 1][0]
                if before in (STATE_BELOW, STATE_ABOVE) and after == -before:
                    events.append({"molecule": i, "direction": "up" if before == STATE_BELOW else "down",
                                   "entry_frame": int(a), "exit_frame": int(runs[r + 1][1]), "status": "complete"})
        # Unresolved: below -> above (or reverse) between consecutive frames: the previous
        # image position plus the unwrapped step lands beyond the far end of the extent.
        s = s_unwrapped[:, i]
        start = s[:-1] if s_image is None else s_image[:-1, i]
        landed = start + np.diff(s)
        prev_state = states[:-1, i]
        up = (prev_state == STATE_BELOW) & (landed > zmax)
        down = (prev_state == STATE_ABOVE) & (landed < zmin)
        for k in np.flatnonzero(up | down):
            events.append({"molecule": i, "direction": "up" if up[k] else "down", "entry_frame": int(k),
                           "exit_frame": int(k + 1), "status": "unresolved_single_interval"})
    return events


def reference_collective_coordinate(s_unwrapped, zmin, zmax, period):
    """Zhu's collective coordinate n(t) for molecules that are always laterally inside.

    n(t) = sum_i [G(s_i(t)) - G(s_i(0))] / L with G(u) the length of the part
    of [lo, u] (in unwrapped coordinates, counting every periodic image of the
    channel) that lies inside the channel: each full passage adds exactly L.
    This is the state-function form of dn = sum (clipped axial displacement) / L
    and is valid when no molecule leaves the lumen laterally.
    """
    L = zmax - zmin
    u = np.asarray(s_unwrapped, dtype=float)
    m = np.floor((u - zmin) / period)
    r = u - zmin - m * period
    g = m * L + np.clip(r, 0, L)
    return (g - g[0]).sum(axis=1) / L


def summarize_events(events):
    """Counts per direction of complete events, unresolved crossings and transits outside the lumen."""
    complete = [e for e in events if e["status"] == "complete"]
    up = sum(e["direction"] == "up" for e in complete)
    return {"up": up, "down": len(complete) - up,
            "unresolved": sum(e["status"] == "unresolved_single_interval" for e in events),
            "outside_lumen": sum(e["status"] == "outside_lumen" for e in events)}


def event_keys(events, offset=0):
    """Hashable (molecule, direction, entry, exit, status) tuples for set comparison."""
    return sorted((int(e["molecule"]) - offset, e["direction"], int(e["entry_frame"]), int(e["exit_frame"]),
                   e["status"]) for e in events)


def crevice_event_keys(result, n_protein, frame_numbers=None):
    """The same tuples from a CREVICE PermeationResult (molecule = atom index - protein atoms).

    CREVICE reports trajectory frame numbers; ``frame_numbers`` maps them back to analysed-frame indices.
    """
    index = None if frame_numbers is None else {f: k for k, f in enumerate(frame_numbers)}
    keys = []
    for e in result.events:
        entry, exit_ = int(e["entry_frame"]), int(e["exit_frame"])
        if index is not None:
            entry, exit_ = index[entry], index[exit_]
        keys.append((int(e["atom_index"]) - n_protein, e["direction"], entry, exit_, e["status"]))
    return sorted(keys)


def dump(path, data):
    """Write JSON with NumPy scalars converted."""
    def convert(value):
        if isinstance(value, dict):
            return {str(k): convert(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [convert(v) for v in value]
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return value
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(convert(data), indent=2) + "\n", encoding="utf-8")
    return Path(path)
