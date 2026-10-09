from __future__ import annotations

import base64

import numpy as np
import pytest

from experiments.gdrive_report.compare import compare_values, plot_traces, scientific_trace


def test_report_comparison_parses_generated_plotly_calls(tmp_path):
    html = tmp_path / "report.html"
    html.write_text('Plotly.newPlot( "random-id", [{"x":[1,2],"y":[3.0,4.0]}], {}, {});')
    assert list(plot_traces(html)) == [[{"x": [1, 2], "y": [3.0, 4.0]}]]


def test_report_comparison_checks_binary_counts_exactly_and_float_rounding():
    def encoded(values):
        return {"dtype": "i4", "bdata": base64.b64encode(np.asarray(values, dtype="i4").tobytes()).decode()}
    compare_values(encoded([1, 2]), encoded([1, 2]), "counts")
    with pytest.raises(AssertionError):
        compare_values(encoded([1, 2]), encoded([1, 3]), "counts")
    compare_values(1., 1. + 1e-15, "quantile")
    with pytest.raises(AssertionError):
        compare_values(1., 1.01, "quantile")


def test_report_comparison_ignores_only_violin_presentation():
    old = {"type": "violin", "name": "BWA", "y": [1.0, 2.0], "jitter": 0.25, "points": "all"}
    new = {"type": "violin", "name": "BWA<br>n=2", "y": [1.0, 2.0],
           "jitter": 0, "points": False, "meanline": {"visible": True}}
    compare_values(scientific_trace(old), scientific_trace(new), "violin")
    new["y"] = [1.0, 3.0]
    with pytest.raises(AssertionError):
        compare_values(scientific_trace(old), scientific_trace(new), "violin")
