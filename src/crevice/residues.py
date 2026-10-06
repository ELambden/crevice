"""Residue-level summaries of region contacts.

The contacts themselves are computed by
:func:`crevice.analysis.annotate_residues` (re-exported here with
:func:`crevice.analysis.nearest_residues`), which measures the gap between
each residue's atom surfaces and a sampled region (a profile, tunnel or
cast). This module turns the resulting :class:`~crevice.models.ResidueContact` records into
tables and label lists. The command-line equivalent is ``crevice residues``.

Roles are geometric: ``"lining"`` means an atom surface within 1 Å of the
region surface, ``"bottleneck"`` additionally means within two samples of the
narrowest profile point. Neither is evidence that a residue controls
transport.
"""

from __future__ import annotations

from typing import Sequence

from .analysis import annotate_residues, nearest_residues
from .models import ResidueContact
from .radii import RadiusSet, radii_option


def residue_contact_table(contacts: Sequence[ResidueContact]) -> list[dict[str, object]]:
    """Convert contacts to plain dictionaries for JSON or CSV output.

    Parameters
    ----------
    contacts : sequence of ResidueContact

    Returns
    -------
    list of dict
        One :meth:`~crevice.models.ResidueContact.to_dict` row per contact, in input order.
    """

    return [contact.to_dict() for contact in contacts]


@radii_option
def lining_residues(
    contacts: Sequence[ResidueContact],
    *,
    include_nearby: bool = False,
    radii: RadiusSet | str | None = None,
) -> tuple[str, ...]:
    """Labels of residues whose role marks them as lining the region.

    Parameters
    ----------
    contacts : sequence of ResidueContact
    include_nearby : bool, default False
        Also include ``"nearby"`` and ``"bottleneck-nearby"`` roles (atom
        surfaces within the annotation cutoff but more than 1 Å from the region
        surface).
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
    tuple of str
        Residue labels with role ``"lining"`` or ``"bottleneck"`` (plus the
        nearby roles if requested), in input order.

    Examples
    --------
    >>> from crevice.models import ResidueContact
    >>> from crevice.residues import lining_residues, bottleneck_residues
    >>> contacts = [ResidueContact("A:LEU10", 0.2, 0.5, 8, 3, "bottleneck"),
    ...             ResidueContact("A:SER14", 0.8, 1.1, 6, 2, "lining"),
    ...             ResidueContact("A:GLU20", 3.0, 3.5, 9, 1, "nearby")]
    >>> lining_residues(contacts)
    ('A:LEU10', 'A:SER14')
    >>> bottleneck_residues(contacts)
    ('A:LEU10',)
    """

    accepted = {"lining", "bottleneck"}
    if include_nearby:
        accepted.update({"nearby", "bottleneck-nearby"})
    return tuple(contact.residue for contact in contacts if contact.role in accepted)


def bottleneck_residues(contacts: Sequence[ResidueContact]) -> tuple[str, ...]:
    """Labels of residues with a bottleneck role.

    Returns every contact whose role contains ``"bottleneck"``, i.e. both
    ``"bottleneck"`` and ``"bottleneck-nearby"``, in input order.

    Parameters
    ----------
    contacts : sequence of ResidueContact

    Returns
    -------
    tuple of str
    """

    return tuple(contact.residue for contact in contacts if "bottleneck" in contact.role)
