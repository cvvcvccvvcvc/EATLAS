from __future__ import annotations

import pandas as pd

from analytics.strategy_report import conservation_bin_detail_table, format_table_dataframe


def test_phyloP_quantiles_are_not_formatted_as_percentages() -> None:
    frame = pd.DataFrame(
        {
            "Background median Q2.5": [0.095],
            "Comparator rate Q2.5": [0.095],
        }
    )

    shown = format_table_dataframe(frame)

    assert shown.loc[0, "Background median Q2.5"] == "0.095"
    assert shown.loc[0, "Comparator rate Q2.5"] == "9.5%"


def test_conservation_bin_table_suppresses_unestimable_inference() -> None:
    bins = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "bin_index": 1,
                "bin_label": "Central phyloP band",
                "bin_range": "-1.30103 to 1.30103",
                "row_count": 10,
                "benign_observed": 6,
                "pathogenic_observed": 0,
                "benign_not_observed": 4,
                "pathogenic_not_observed": 0,
                "odds_ratio": float("nan"),
                "ci_low": 0.1,
                "ci_high": 10.0,
                "fisher_p": 1.0,
                "fisher_q": 1.0,
            }
        ]
    )

    table = conservation_bin_detail_table(bins)

    assert table.loc[0, "95% CI"] == ""
    assert table.loc[0, "Fisher p"] == ""
    assert table.loc[0, "Status"] == "Not estimable"
