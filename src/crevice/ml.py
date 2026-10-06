"""Simple, interpretable feature tables, clustering and state summaries.

These utilities turn profiles and residue contacts into named numeric
features, standardise them, and apply small deterministic methods: k-means
clustering, nearest-to-centroid representatives, PCA projection and
transition counting between cluster labels. They are dependency-free except
that PCA uses scikit-learn when it is installed.

Standardisation
    Every method that compares feature rows first takes the union of feature
    names (missing values become 0.0) and converts each column to z-scores
    (subtract the mean, divide by the population standard deviation; a
    constant column is divided by 1).

Limitations
    k-means starts from the first ``k`` rows, so results depend on input
    order and can reach a local optimum. Cluster labels are arbitrary
    integers, and clusters of geometric features are descriptive groupings,
    not thermodynamic states.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

from .models import PoreProfile, ResidueContact


def profile_features(
    profile: PoreProfile,
    contacts: Sequence[ResidueContact] | None = None,
) -> dict[str, float]:
    """Named numeric features of one profile (and optionally its contacts).

    Parameters
    ----------
    profile : PoreProfile
    contacts : sequence of ResidueContact, optional
        Adds contact-composition features when given.

    Returns
    -------
    dict of str to float
        Profile features: ``min_radius``, ``max_radius``, ``mean_radius``
        (Å), ``length`` (Å), ``volume_estimate`` (Å³), ``bottleneck_position``
        (``t`` of the narrowest point) and ``radius_variance`` (population
        variance, Å²). With contacts also: ``contact_count``,
        ``hydrophobic_fraction``, ``polar_fraction``, ``charged_fraction``
        (positive plus negative classes, so HIS counts as charged),
        ``bottleneck_contact_count`` and ``mean_residue_influence``. Fractions
        and the mean are over the number of contacts (1 if there are none).

    Raises
    ------
    ValueError
        If the profile has no points.
    """

    features = {
        "min_radius": profile.min_radius,
        "max_radius": profile.max_radius,
        "mean_radius": profile.mean_radius,
        "length": profile.length,
        "volume_estimate": profile.volume_estimate,
        "bottleneck_position": profile.bottleneck.t,
        "radius_variance": _variance([point.radius for point in profile.points]),
    }
    if contacts is not None:
        total = max(len(contacts), 1)
        features.update(
            {
                "contact_count": float(len(contacts)),
                "hydrophobic_fraction": _property_fraction(contacts, "hydrophobic", total),
                "polar_fraction": _property_fraction(contacts, "polar", total),
                "charged_fraction": (
                    _property_fraction(contacts, "positive", total)
                    + _property_fraction(contacts, "negative", total)
                ),
                "bottleneck_contact_count": float(sum("bottleneck" in contact.role for contact in contacts)),
                "mean_residue_influence": sum(contact.influence_score for contact in contacts) / total,
            }
        )
    return features


def representative_profile(
    profiles: Sequence[PoreProfile],
    contacts: Sequence[Sequence[ResidueContact] | None] | None = None,
) -> tuple[int, PoreProfile]:
    """Profile closest to the centroid of the set in standardised feature space.

    Parameters
    ----------
    profiles : sequence of PoreProfile
    contacts : sequence of (sequence of ResidueContact or None), optional
        Per-profile contacts, aligned with ``profiles``.

    Returns
    -------
    index : int
        Position of the representative (the first one on a tie).
    profile : PoreProfile

    Raises
    ------
    ValueError
        If ``profiles`` is empty.
    """

    if not profiles:
        raise ValueError("At least one profile is required")
    contact_sets: Sequence[Sequence[ResidueContact] | None]
    contact_sets = contacts if contacts is not None else [None] * len(profiles)
    rows = [profile_features(profile, contact) for profile, contact in zip(profiles, contact_sets)]
    keys = sorted({key for row in rows for key in row})
    matrix = [[row.get(key, 0.0) for key in keys] for row in rows]
    normalized = _normalize_columns(matrix)
    center = [sum(row[col] for row in normalized) / len(normalized) for col in range(len(keys))]
    distances = [_euclidean(row, center) for row in normalized]
    index = min(range(len(distances)), key=distances.__getitem__)
    return index, profiles[index]


def cluster_profiles(
    profiles: Sequence[PoreProfile],
    *,
    k: int = 2,
    iterations: int = 50,
) -> list[int]:
    """Cluster profiles by k-means on standardised profile features.

    Contact features are not used. Initial centres are the first ``k`` rows;
    points are assigned to the nearest centre (lowest label on ties) and centres
    recomputed until no label changes or ``iterations`` is reached. An empty
    cluster keeps its previous centre.

    Parameters
    ----------
    profiles : sequence of PoreProfile
    k : int, default 2
        Number of clusters (reduced to the number of profiles if larger).
    iterations : int, default 50
        Maximum iterations.

    Returns
    -------
    list of int
        Label (0 to ``k - 1``) for each profile; empty for no profiles.

    Raises
    ------
    ValueError
        If ``k < 1``.
    """

    if k < 1:
        raise ValueError("k must be at least 1")
    if not profiles:
        return []
    rows = [profile_features(profile) for profile in profiles]
    keys = sorted({key for row in rows for key in row})
    matrix = _normalize_columns([[row.get(key, 0.0) for key in keys] for row in rows])
    k = min(k, len(matrix))
    centers = [matrix[index][:] for index in range(k)]
    labels = [0] * len(matrix)
    for _ in range(iterations):
        changed = False
        for row_index, row in enumerate(matrix):
            label = min(range(k), key=lambda center_index: _euclidean(row, centers[center_index]))
            if labels[row_index] != label:
                labels[row_index] = label
                changed = True
        new_centers: list[list[float]] = []
        for center_index in range(k):
            members = [row for row, label in zip(matrix, labels) if label == center_index]
            if not members:
                new_centers.append(centers[center_index])
                continue
            new_centers.append(
                [sum(row[col] for row in members) / len(members) for col in range(len(keys))]
            )
        centers = new_centers
        if not changed:
            break
    return labels


def residue_importance(contacts: Sequence[ResidueContact]) -> list[dict[str, object]]:
    """Contacts as dictionaries, sorted by influence score (highest first), then
    by smallest ``min_distance``.

    The influence score is a geometric proximity heuristic (see
    :class:`~crevice.models.ResidueContact`), not a functional importance.

    Parameters
    ----------
    contacts : sequence of ResidueContact
        From :func:`crevice.analysis.annotate_residues`.

    Returns
    -------
    list of dict
        :meth:`~crevice.models.ResidueContact.to_dict` of each contact, in rank
        order.
    """

    return [
        contact.to_dict()
        for contact in sorted(contacts, key=lambda item: (-item.influence_score, item.min_distance))
    ]


def _property_fraction(contacts: Sequence[ResidueContact], prop: str, total: int) -> float:
    """Fraction of contacts whose residue has class ``prop``."""
    return sum(prop in contact.properties for contact in contacts) / total


def _variance(values: Sequence[float]) -> float:
    """Population variance; 0.0 for an empty sequence."""
    if not values:
        return 0.0
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / len(values)


def _normalize_columns(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    """Column-wise z-scores (population standard deviation; constant columns divided by 1)."""
    if not matrix:
        return []
    cols = len(matrix[0])
    means = [sum(row[col] for row in matrix) / len(matrix) for col in range(cols)]
    variances = [
        sum((row[col] - means[col]) ** 2 for row in matrix) / len(matrix)
        for col in range(cols)
    ]
    scales = [math.sqrt(value) or 1.0 for value in variances]
    return [
        [(row[col] - means[col]) / scales[col] for col in range(cols)]
        for row in matrix
    ]


def _euclidean(left: Sequence[float], right: Sequence[float]) -> float:
    """Euclidean distance between two equal-length vectors."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def extract_features(
    results: Sequence[object] | object,
    *,
    include_contacts: bool = True,
) -> list[dict[str, float]]:
    """Feature rows from profiles, frame results or a trajectory result.

    Parameters
    ----------
    results : object or sequence
        A :class:`~crevice.models.TrajectoryAnalysis` (anything with ``frames``), a sequence of
        items, or a single item. Each item may be a :class:`~crevice.models.PoreProfile`, an
        object with a non-empty ``features`` dict (used as is), an object with
        a ``profile`` (features computed with :func:`profile_features`), or a
        mapping (its ``int``/``float`` values are kept). Other items are
        skipped.
    include_contacts : bool, default True
        Include contact features from an item's ``residue_contacts``.

    Returns
    -------
    list of dict of str to float
        One row per usable item.
    """

    if hasattr(results, "frames"):
        iterable = list(getattr(results, "frames"))
    elif isinstance(results, Sequence) and not isinstance(results, (str, bytes)):
        iterable = list(results)
    else:
        iterable = [results]

    rows: list[dict[str, float]] = []
    for item in iterable:
        if isinstance(item, PoreProfile):
            rows.append(profile_features(item))
            continue
        profile = getattr(item, "profile", None)
        contacts = getattr(item, "residue_contacts", None) if include_contacts else None
        existing = getattr(item, "features", None)
        if existing:
            rows.append({key: float(value) for key, value in existing.items()})
        elif profile is not None:
            rows.append(profile_features(profile, contacts))
        elif isinstance(item, Mapping):
            rows.append({str(key): float(value) for key, value in item.items() if isinstance(value, (int, float))})
    return rows


def reduce_dimensions(
    features: Sequence[Mapping[str, float]],
    *,
    method: str = "pca",
    components: int = 2,
) -> list[dict[str, float]]:
    """Project standardised feature rows to a few coordinates.

    Parameters
    ----------
    features : sequence of mapping of str to float
    method : {"pca", "identity"}, default "pca"
        ``"pca"`` uses :class:`sklearn.decomposition.PCA` when scikit-learn is
        installed; otherwise it silently falls back to ``"identity"``.
        ``"identity"`` keeps the first standardised columns in alphabetical
        feature-name order.
    components : int, default 2
        Output dimensions; missing dimensions are padded with 0.0.

    Returns
    -------
    list of dict
        Rows with keys ``dim_1``, ``dim_2``, ...

    Raises
    ------
    ValueError
        If ``components < 1`` or the method is unknown.

    Examples
    --------
    >>> from crevice.ml import reduce_dimensions
    >>> reduce_dimensions([{"a": 1.0, "b": 5.0}, {"a": 3.0, "b": 5.0}], method="identity")
    [{'dim_1': -1.0, 'dim_2': 0.0}, {'dim_1': 1.0, 'dim_2': 0.0}]
    """

    if components < 1:
        raise ValueError("components must be at least 1")
    if not features:
        return []
    method = method.lower()
    keys = sorted({key for row in features for key in row})
    matrix = _normalize_columns([[float(row.get(key, 0.0)) for key in keys] for row in features])
    if method == "identity":
        projected = [row[:components] + [0.0] * max(0, components - len(row)) for row in matrix]
    elif method == "pca":
        projected = _pca_or_fallback(matrix, components)
    else:
        raise ValueError("method must be identity or pca in the dependency-free core")
    return [
        {f"dim_{index + 1}": value for index, value in enumerate(row)}
        for row in projected
    ]


def cluster_conformations(
    features: Sequence[Mapping[str, float]],
    *,
    k: int = 2,
    iterations: int = 50,
    method: str = "kmeans",
) -> list[int]:
    """Cluster arbitrary feature rows with deterministic k-means.

    Same algorithm as :func:`cluster_profiles`, applied to the standardised
    union of feature names.

    Parameters
    ----------
    features : sequence of mapping of str to float
    k : int, default 2
    iterations : int, default 50
    method : {"kmeans"}, default "kmeans"

    Returns
    -------
    list of int

    Raises
    ------
    ValueError
        If ``k < 1`` or the method is not ``"kmeans"``.

    Examples
    --------
    >>> from crevice.ml import cluster_conformations
    >>> rows = [{"r": 1.0}, {"r": 1.1}, {"r": 4.0}, {"r": 4.2}]
    >>> cluster_conformations(rows, k=2)
    [0, 0, 1, 1]
    """

    if method != "kmeans":
        raise ValueError("method must be kmeans in the dependency-free core")
    if k < 1:
        raise ValueError("k must be at least 1")
    if not features:
        return []
    keys = sorted({key for row in features for key in row})
    matrix = _normalize_columns([[float(row.get(key, 0.0)) for key in keys] for row in features])
    k = min(k, len(matrix))
    centers = [matrix[index][:] for index in range(k)]
    labels = [0] * len(matrix)
    for _ in range(iterations):
        changed = False
        for row_index, row in enumerate(matrix):
            label = min(range(k), key=lambda center_index: _euclidean(row, centers[center_index]))
            if labels[row_index] != label:
                labels[row_index] = label
                changed = True
        new_centers = []
        for center_index in range(k):
            members = [row for row, label in zip(matrix, labels) if label == center_index]
            if not members:
                new_centers.append(centers[center_index])
            else:
                new_centers.append([sum(row[col] for row in members) / len(members) for col in range(len(keys))])
        centers = new_centers
        if not changed:
            break
    return labels


def representative_frames(
    features: Sequence[Mapping[str, float]],
    labels: Sequence[int] | None = None,
) -> dict[int, int]:
    """Row closest to each cluster's centroid in standardised feature space.

    Parameters
    ----------
    features : sequence of mapping of str to float
    labels : sequence of int, optional
        Cluster label per row; all rows form one cluster when omitted.

    Returns
    -------
    dict of int to int
        Cluster label to row index (the first row on a tie).

    Raises
    ------
    ValueError
        If ``labels`` and ``features`` differ in length.
    """

    if not features:
        return {}
    labels = list(labels) if labels is not None else [0] * len(features)
    if len(labels) != len(features):
        raise ValueError("labels must have the same length as features")
    keys = sorted({key for row in features for key in row})
    matrix = _normalize_columns([[float(row.get(key, 0.0)) for key in keys] for row in features])
    representatives: dict[int, int] = {}
    for label in sorted(set(labels)):
        indices = [index for index, value in enumerate(labels) if value == label]
        center = [sum(matrix[index][col] for index in indices) / len(indices) for col in range(len(keys))]
        representatives[label] = min(indices, key=lambda index: _euclidean(matrix[index], center))
    return representatives


def transition_summary(labels: Sequence[int], times: Sequence[float] | None = None) -> dict[str, object]:
    """Count label changes between consecutive frames.

    Parameters
    ----------
    labels : sequence of int
        One state label per frame, in time order.
    times : sequence of float, optional
        Frame times, aligned with ``labels``.

    Returns
    -------
    dict
        ``state_counts`` (frames per label), ``transition_counts`` (keys
        ``"a->b"``) and ``events``: one per change, with ``index`` of the frame
        after the change, its ``time`` (or ``None``) and ``from``/``to`` labels.

    Notes
    -----
    Counts depend on the frame spacing; changes between saved frames are not
    seen, and no rates or lifetimes are fitted.

    Examples
    --------
    >>> from crevice.ml import transition_summary
    >>> summary = transition_summary([0, 0, 1, 1, 0])
    >>> summary["state_counts"], summary["transition_counts"]
    ({0: 3, 1: 2}, {'0->1': 1, '1->0': 1})
    """

    transitions: dict[str, int] = {}
    dwell: dict[int, int] = {}
    for label in labels:
        dwell[label] = dwell.get(label, 0) + 1
    events = []
    for index, (left, right) in enumerate(zip(labels, labels[1:]), start=1):
        if left == right:
            continue
        key = f"{left}->{right}"
        transitions[key] = transitions.get(key, 0) + 1
        events.append({"index": index, "time": times[index] if times is not None else None, "from": left, "to": right})
    return {"state_counts": dwell, "transition_counts": transitions, "events": events}


def _pca_or_fallback(matrix: Sequence[Sequence[float]], components: int) -> list[list[float]]:
    """PCA via scikit-learn, or the first standardised columns if it is not installed."""
    try:
        from sklearn.decomposition import PCA  # type: ignore[import-not-found]
    except ImportError:
        return [list(row[:components]) + [0.0] * max(0, components - len(row)) for row in matrix]
    projected = PCA(n_components=min(components, len(matrix[0]), len(matrix))).fit_transform(matrix)
    rows = [list(map(float, row)) for row in projected]
    for row in rows:
        row.extend([0.0] * max(0, components - len(row)))
    return [row[:components] for row in rows]

