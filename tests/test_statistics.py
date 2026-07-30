from __future__ import annotations

import pytest

from analytics.analyses.statistics import fisher_exact_two_sided


@pytest.mark.parametrize(
    ("cells", "expected"),
    [
        ((1, 9, 11, 3), 0.0027594561852200836),
        ((0, 5, 4, 1), 0.04761904761904762),
        ((5, 5, 5, 5), 1.0),
    ],
)
def test_fisher_exact_two_sided_matches_reference_values(
    cells: tuple[int, int, int, int],
    expected: float,
) -> None:
    assert fisher_exact_two_sided(*cells) == pytest.approx(expected)
