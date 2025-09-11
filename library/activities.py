from __future__ import annotations


def get_activities(limit: int) -> list[dict[str, int]]:
    """Return dummy activity identifiers.

    Parameters
    ----------
    limit
        Number of activity rows to generate. Must be non-negative.

    Returns
    -------
    list of dict
        Sequence of dictionaries with an ``activity_id`` field.

    Raises
    ------
    ValueError
        If ``limit`` is negative.
    """
    if limit < 0:
        raise ValueError("limit must be non-negative")
    return [{"activity_id": i} for i in range(limit)]
