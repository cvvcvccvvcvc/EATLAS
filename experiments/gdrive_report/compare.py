"""Compare scientific outputs of two completed reports, excluding provenance/UI."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.analyses.variant_summary import _frame_from_payload
from analytics.io.artifacts import write_json_atomic


FLOAT_TOLERANCE = 1e-12
TABLES = (
    "candidate_variants.phyloP100way.distributions.tsv.gz",
    "candidate_variants.phyloP100way.histograms.tsv.gz",
    "clinvar_universe.snv_indel.tsv.gz",
    "clinvar_universe.snv_indel.vep.tsv.gz",
    "clinvar_universe.snv_indel.conservation.tsv.gz",
    "clinvar_observed_memberships.tsv.gz",
    "continuous_firth/results.tsv.gz",
    "continuous_firth/distributions.tsv.gz",
    "pathogenic_clinvar_hits.tsv.gz",
)


def sorted_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.reset_index(drop=True)
    return frame.sort_values(list(frame.columns), kind="stable").reset_index(drop=True)


def compare_values(left, right, label: str) -> None:
    if isinstance(left, dict) and isinstance(right, dict):
        # Plotly's random trace IDs are presentation state, not measurements.
        keys = set(left) - {"uid"}
        assert keys == set(right) - {"uid"}, f"Different fields: {label}"
        if "bdata" in keys and "dtype" in keys:
            assert left["dtype"] == right["dtype"], f"Different dtype: {label}"
            dtype = np.dtype(left["dtype"])
            a = np.frombuffer(base64.b64decode(left["bdata"]), dtype=dtype)
            b = np.frombuffer(base64.b64decode(right["bdata"]), dtype=dtype)
            if np.issubdtype(dtype, np.integer):
                np.testing.assert_array_equal(a, b, err_msg=label)
            else:
                np.testing.assert_allclose(a, b, rtol=FLOAT_TOLERANCE, atol=FLOAT_TOLERANCE, err_msg=label)
            keys -= {"bdata", "dtype"}
        for key in sorted(keys):
            compare_values(left[key], right[key], f"{label}.{key}")
    elif isinstance(left, list) and isinstance(right, list):
        assert len(left) == len(right), f"Different lengths: {label}"
        for index, (a, b) in enumerate(zip(left, right)):
            compare_values(a, b, f"{label}[{index}]")
    elif isinstance(left, float) or isinstance(right, float):
        assert np.isclose(left, right, rtol=FLOAT_TOLERANCE, atol=FLOAT_TOLERANCE, equal_nan=True), label
    else:
        assert left == right, f"Different value: {label}"


def plot_traces(html_path: Path):
    html = html_path.read_text()
    decoder = json.JSONDecoder()
    for match in re.finditer(r"Plotly\.newPlot\(\s*(?=[\"'])", html):
        # Generated reports use a JSON string ID followed by JSON trace data.
        start = match.end()
        _, end = decoder.raw_decode(html, start)
        start = end + 1
        while html[start].isspace():
            start += 1
        traces, _ = decoder.raw_decode(html, start)
        yield traces


def scientific_trace(trace: dict) -> dict:
    if trace.get("type") != "violin":
        return trace
    # These fields change the violin's presentation, not its plotted values.
    result = {key: value for key, value in trace.items()
              if key not in {"hoveron", "jitter", "meanline", "points"}}
    if isinstance(result.get("name"), str):
        result["name"] = re.sub(r"<br>n=[\d,]+$", "", result["name"])
    return result


def compare_reports(baseline: Path, candidate: Path, baseline_html: Path, candidate_html: Path) -> dict:
    manifests = [json.loads((directory / "manifest.json").read_text()) for directory in (baseline, candidate)]
    assert {source["source_id"] for source in manifests[0]["sources"]} == {
        source["source_id"] for source in manifests[1]["sources"]
    }, "Different scientific source cohorts"
    with gzip.open(baseline / "derived/variant_summary.json.gz", "rt") as handle:
        old = json.load(handle)["summary"]
    with gzip.open(candidate / "derived/variant_summary.json.gz", "rt") as handle:
        new = json.load(handle)["summary"]
    compare_values(
        {k: v for k, v in old.items() if k != "frames"},
        {k: v for k, v in new.items() if k != "frames"}, "summary",
    )
    assert set(old["frames"]) == set(new["frames"]), "Different scientific summary tables"
    for name, payload in old["frames"].items():
        pd.testing.assert_frame_equal(
            sorted_frame(_frame_from_payload(payload)),
            sorted_frame(_frame_from_payload(new["frames"][name])),
            rtol=FLOAT_TOLERANCE, atol=FLOAT_TOLERANCE, obj=name,
        )
    for relative in TABLES:
        frames = [pd.read_csv(directory / "derived" / relative, sep="\t", keep_default_na=False)
                  for directory in (baseline, candidate)]
        pd.testing.assert_frame_equal(
            sorted_frame(frames[0]), sorted_frame(frames[1]),
            rtol=FLOAT_TOLERANCE, atol=FLOAT_TOLERANCE, obj=relative,
        )
    from itertools import zip_longest

    plot_count = 0
    for index, (old_plot, new_plot) in enumerate(zip_longest(plot_traces(baseline_html), plot_traces(candidate_html))):
        assert old_plot is not None and new_plot is not None, "Different number of scientific figures"
        compare_values(
            [scientific_trace(trace) for trace in old_plot],
            [scientific_trace(trace) for trace in new_plot],
            f"report.scientific_traces[{index}]",
        )
        plot_count += 1
    assert plot_count > 0, "No report figures parsed; comparison is incomplete"
    return {"status": "identical_scientific_results", "summary_tables": len(old["frames"]),
            "derived_tables": len(TABLES), "plot_count": plot_count,
            "float_tolerance": FLOAT_TOLERANCE,
            "baseline": str(baseline), "candidate": str(candidate)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-analysis-dir", required=True, type=Path)
    parser.add_argument("--baseline-report-name", required=True)
    parser.add_argument("--analytics-root", required=True, type=Path)
    parser.add_argument("--candidate-report-name", required=True)
    parser.add_argument("--result", required=True, type=Path)
    args = parser.parse_args()
    candidates = list(args.analytics_root.glob(f"analyses/*/reports/{args.candidate_report_name}.html"))
    if len(candidates) != 1:
        raise ValueError(f"Expected one completed candidate report, found {len(candidates)}")
    candidate_html = candidates[0]
    try:
        result = compare_reports(
            args.baseline_analysis_dir, candidate_html.parent.parent,
            args.baseline_analysis_dir / "reports" / f"{args.baseline_report_name}.html", candidate_html,
        )
    except Exception as error:
        write_json_atomic(args.result, {"status": "failed", "error_type": type(error).__name__,
                                       "error": str(error)[:1500]})
        raise
    write_json_atomic(args.result, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
