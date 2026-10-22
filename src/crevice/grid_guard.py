"""Warn when a grid is too coarse for the narrowest clearance it has to resolve.

Grid-based measurements (casts, cavity-trajectory regions, path searches) lose
narrow openings and volume as the grid coarsens. On 1GRM the rolling-probe cast
volume falls from 289.9 to 282.8, 257.0 and 238.4 Å³ at 0.25, 0.35, 0.5 and
0.7 Å (-2.4 / -11.3 / -17.8 % against 0.25 Å), and synthetic walls lost their
narrow opening at coarse grids (CREVICE 0.3 development runs). These
are software/benchmark observations, not a validated error model.

:func:`grid_resolution_check` flags a run whose narrowest clearance is below
about two grid spacings and recommends a spacing of at most 0.35 Å for narrow
channels. It never changes a setting or a result; commands record the check in
their manifest (``grid_resolution``) and print a warning.

What counts as the narrowest clearance
    * a resolved channel profile: its minimum radius (``publish``);
    * a path search: the narrowest bottleneck radius of the reported routes.

Only measured clearances are checked. A rolling-probe cast or per-frame region
(``cavities --cast``, ``trajectory --cavities``, ``trajectory --region``)
reports no narrowest width, so its check is recorded as ``not_assessed`` and
no warning is printed.
"""

from __future__ import annotations

import math

#: Clearances below this many grid spacings trigger the warning.
GRID_GUARD_FACTOR = 2.0
#: Spacing recommended for narrow channels (Å).
RECOMMENDED_SPACING_A = 0.35
EVIDENCE = ("1GRM rolling-probe cast volume 289.9 / 282.8 / 257.0 / 238.4 A^3 at grid 0.25 / 0.35 / 0.5 / 0.7 A; "
            "synthetic walls lost narrow openings at coarse grids (CREVICE 0.3 development runs)")


def grid_resolution_check(spacing_A, narrowest_clearance_A, *, clearance_source: str, context: str = "") -> dict:
    """Compare the grid spacing with the narrowest clearance the run must resolve.

    Parameters
    ----------
    spacing_A : float
        Grid spacing used (Å).
    narrowest_clearance_A : float or None
        Narrowest clearance (Å) of the route or cast; ``None`` when unknown.
    clearance_source : str
        What the clearance is (recorded verbatim).
    context : str, optional
        What was measured, for the message (for example ``"cast"``).

    Returns
    -------
    dict
        ``status`` (``"warning"``, ``"ok"`` or ``"not_assessed"``),
        ``spacing_A``, ``narrowest_clearance_A``, ``clearance_source``,
        ``threshold_A`` (``2 * spacing``), ``recommended_spacing_A``,
        ``message`` (``None`` unless warning), ``rule`` and ``evidence``.

    Examples
    --------
    >>> from crevice.grid_guard import grid_resolution_check
    >>> grid_resolution_check(0.7, 1.33, clearance_source="profile minimum radius")["status"]
    'warning'
    >>> grid_resolution_check(0.35, 1.33, clearance_source="profile minimum radius")["status"]
    'ok'
    """
    spacing = float(spacing_A)
    threshold = GRID_GUARD_FACTOR * spacing
    record = {"status": "not_assessed", "spacing_A": spacing, "narrowest_clearance_A": None,
              "clearance_source": clearance_source, "threshold_A": threshold,
              "recommended_spacing_A": RECOMMENDED_SPACING_A, "message": None,
              "rule": f"warn when the narrowest clearance is below {GRID_GUARD_FACTOR:g} grid spacings; "
                      "settings and results are never changed",
              "evidence": EVIDENCE}
    if narrowest_clearance_A is None or not math.isfinite(float(narrowest_clearance_A)):
        return record
    clearance = float(narrowest_clearance_A)
    record["narrowest_clearance_A"] = clearance
    if clearance < threshold:
        what = f" of the {context}" if context else ""
        suggestion = (f"use a spacing of at most {RECOMMENDED_SPACING_A:g} A for narrow channels" if spacing > RECOMMENDED_SPACING_A
                      else "use a finer spacing")
        record.update(status="warning", message=(
            f"the narrowest clearance{what} ({clearance:.2f} A, {clearance_source}) is below two grid spacings "
            f"({threshold:.2f} A at {spacing:g} A); narrow openings and volume may be lost: {suggestion}"))
    else:
        record["status"] = "ok"
    return record
