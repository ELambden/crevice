"""Versioned envelope for CREVICE JSON outputs.

Every JSON result written by CREVICE can be wrapped in the same envelope so
that readers can check what kind of result a file holds and which schema
version produced it. :data:`SCHEMA_VERSION` is the current version string.
"""

from __future__ import annotations

from typing import Any, Mapping

#: Version string embedded in every result envelope schema id.
SCHEMA_VERSION = "0.1"


def result_envelope(
    kind: str,
    data: Mapping[str, Any] | list[Any],
    *,
    parameters: Mapping[str, Any] | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Wrap a result payload with stable schema metadata.

    Parameters
    ----------
    kind : str
        Result kind, for example ``"profile"``; becomes part of the schema id.
    data : mapping or list
        The payload, stored unchanged.
    parameters : mapping, optional
        Settings used to compute the result (copied into a new dict).
    metadata : mapping, optional
        Provenance or other context (copied into a new dict).

    Returns
    -------
    dict
        ``schema`` (``"crevice.<kind>.v<SCHEMA_VERSION>"``), ``kind``,
        ``parameters``, ``metadata`` and ``data``.

    Examples
    --------
    >>> from crevice.schemas import result_envelope
    >>> result_envelope("profile", {"min_radius": 1.2})["schema"]
    'crevice.profile.v0.1'
    """

    return {
        "schema": f"crevice.{kind}.v{SCHEMA_VERSION}",
        "kind": kind,
        "parameters": dict(parameters or {}),
        "metadata": dict(metadata or {}),
        "data": data,
    }
