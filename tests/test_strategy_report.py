from __future__ import annotations

import pandas as pd

from analytics.strategy_report import conservation_selector_view, dataframe_records, format_table_dataframe


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


def test_conservation_selector_serializes_sparse_results_and_has_all_controls() -> None:
    primary = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "variant_type": "snv",
                "consequence": "missense",
                "odds_ratio_mh": float("inf"),
                "ci_low": float("nan"),
                "ci_high": float("nan"),
                "cmh_p": float("nan"),
                "cmh_q": float("nan"),
                "usable_rows": 10,
                "status": "not_estimable",
                "reason": "Sparse data",
            }
        ]
    )
    detail = pd.DataFrame()

    html = conservation_selector_view(
        view_id="fixed-test",
        strategies=["s1"],
        primary=primary,
        detail=detail,
        mode="fixed",
    )

    assert 'data-role="strategy"' in html
    assert 'data-role="variant-type"' in html
    assert 'data-role="consequence"' in html
    assert "Missense" in html
    assert "Infinity" not in html
    assert "NaN" not in html


def test_dataframe_records_replaces_nonfinite_values() -> None:
    records = dataframe_records(pd.DataFrame({"value": [1.0, float("inf"), float("nan")]}))
    assert records == [{"value": 1.0}, {"value": "inf"}, {"value": None}]
