#!/usr/bin/env python3
"""Evaluate baseline versus GAPH feature ablations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .io import write_tsv


ID_COLUMNS = {
    "variant_id",
    "gene_id",
    "genomic_accession",
    "genomic_start1",
    "ref",
    "alt",
    "target_start0",
    "target_end0",
    "event_type",
    "target_feature_types",
    "strategy",
}
POSITIVE_LABELS = {"1", "true", "pathogenic", "likely_pathogenic", "p/lp", "lp/p", "p", "lp"}
NEGATIVE_LABELS = {"0", "false", "benign", "likely_benign", "b/lb", "lb/b", "b", "lb"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--baseline-features", default="", help="Comma-separated baseline columns.")
    parser.add_argument("--gaph-prefix", default="gaph_")
    parser.add_argument("--group-column", help="Use group-aware CV on this column, for example gene_id.")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--model", choices=["logistic", "linear_svm"], default="logistic")
    parser.add_argument("--bootstrap-iterations", type=int, default=500)
    parser.add_argument("--random-seed", type=int, default=13)
    return parser.parse_args()


def parse_label(value: object) -> int | None:
    text = str(value).strip().lower().replace(" ", "_")
    if text in POSITIVE_LABELS:
        return 1
    if text in NEGATIVE_LABELS:
        return 0
    return None


def parse_feature_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def numeric_feature_frame(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for col in columns:
        if col not in df.columns:
            raise ValueError(f"Requested feature column not found: {col}")
        out[col] = pd.to_numeric(df[col], errors="coerce")
    return out


def infer_gaph_features(df: pd.DataFrame, prefix: str) -> list[str]:
    return [
        col
        for col in df.columns
        if col.startswith(prefix) and pd.to_numeric(df[col], errors="coerce").notna().any()
    ]


def infer_baseline_features(df: pd.DataFrame, label_column: str, gaph_features: list[str]) -> list[str]:
    excluded = set(gaph_features) | ID_COLUMNS | {label_column}
    out = []
    for col in df.columns:
        if col in excluded:
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        if numeric.notna().any():
            out.append(col)
    return out


def split_indices(df: pd.DataFrame, y: np.ndarray, group_column: str | None, folds: int, seed: int):
    from sklearn.model_selection import GroupKFold, StratifiedGroupKFold, StratifiedKFold

    class_min = int(np.bincount(y).min())
    folds = max(2, min(folds, class_min))
    if group_column:
        if group_column not in df.columns:
            raise ValueError(f"Group column not found: {group_column}")
        groups = df[group_column].astype(str).to_numpy()
        unique_groups = len(set(groups))
        if unique_groups < 2:
            raise ValueError(f"Group column {group_column} has fewer than two groups")
        folds = min(folds, unique_groups)
        try:
            return list(StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed).split(df, y, groups))
        except ValueError:
            return list(GroupKFold(n_splits=folds).split(df, y, groups))
    return list(StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed).split(df, y))


def split_assignment_rows(df: pd.DataFrame, splits: list[tuple[np.ndarray, np.ndarray]]) -> list[dict[str, object]]:
    rows = []
    assignment = {}
    for fold, (_train_idx, test_idx) in enumerate(splits, start=1):
        for row_index in test_idx:
            assignment[int(row_index)] = fold
    for row_index in range(len(df)):
        rows.append(
            {
                "variant_id": df.iloc[row_index].get("variant_id", ""),
                "gene_id": df.iloc[row_index].get("gene_id", ""),
                "strategy": df.iloc[row_index].get("strategy", ""),
                "fold": assignment.get(row_index, ""),
            }
        )
    return rows


def make_model(model_name: str):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import LinearSVC

    if model_name == "logistic":
        estimator = LogisticRegression(max_iter=2000, class_weight="balanced")
    else:
        estimator = LinearSVC(class_weight="balanced", dual="auto", max_iter=5000)
    return Pipeline(
        steps=[
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", estimator),
        ]
    )


def model_scores(model, x_test: pd.DataFrame) -> np.ndarray:
    estimator = model.named_steps["model"]
    if hasattr(estimator, "predict_proba"):
        return model.predict_proba(x_test)[:, 1]
    return model.decision_function(x_test)


def metric_value(y_true: np.ndarray, scores: np.ndarray, metric: str) -> float:
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    if metric == "auroc":
        if len(set(y_true)) != 2:
            return float("nan")
        return float(roc_auc_score(y_true, scores))
    if metric == "auprc":
        return float(average_precision_score(y_true, scores))
    if metric == "brier":
        if len(scores) == 0 or float(np.nanmin(scores)) < 0.0 or float(np.nanmax(scores)) > 1.0:
            return float("nan")
        return float(brier_score_loss(y_true, scores))
    raise ValueError(f"Unknown metric: {metric}")


def bootstrap_ci(
    pred_df: pd.DataFrame,
    metric: str,
    group_column: str | None,
    seed: int,
    iterations: int,
) -> tuple[float, float]:
    if pred_df.empty or iterations <= 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    values = []
    if group_column and group_column in pred_df.columns and pred_df[group_column].astype(str).nunique() > 1:
        groups = sorted(pred_df[group_column].astype(str).unique())
        grouped = {group: pred_df[pred_df[group_column].astype(str) == group] for group in groups}
        for _ in range(iterations):
            sampled = rng.choice(groups, size=len(groups), replace=True)
            boot = pd.concat([grouped[group] for group in sampled], ignore_index=True)
            if boot["label"].nunique() < 2:
                continue
            values.append(metric_value(boot["label"].to_numpy(), boot["score"].to_numpy(), metric))
    else:
        n = len(pred_df)
        for _ in range(iterations):
            idx = rng.integers(0, n, size=n)
            boot = pred_df.iloc[idx]
            if boot["label"].nunique() < 2:
                continue
            values.append(metric_value(boot["label"].to_numpy(), boot["score"].to_numpy(), metric))
    values = [value for value in values if not np.isnan(value)]
    if not values:
        return float("nan"), float("nan")
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def summarize_predictions(
    pred_df: pd.DataFrame,
    fold_metrics: list[dict[str, float]],
    args: argparse.Namespace,
) -> dict[str, object]:
    y_true = pred_df["label"].to_numpy()
    scores = pred_df["score"].to_numpy()
    auroc_ci = bootstrap_ci(pred_df, "auroc", args.group_column, args.random_seed, args.bootstrap_iterations)
    auprc_ci = bootstrap_ci(pred_df, "auprc", args.group_column, args.random_seed, args.bootstrap_iterations)
    brier_ci = bootstrap_ci(pred_df, "brier", args.group_column, args.random_seed, args.bootstrap_iterations)
    return {
        "auroc_oof": metric_value(y_true, scores, "auroc"),
        "auroc_ci_low": auroc_ci[0],
        "auroc_ci_high": auroc_ci[1],
        "auroc_fold_mean": float(np.nanmean([item["auroc"] for item in fold_metrics])),
        "auroc_fold_std": float(np.nanstd([item["auroc"] for item in fold_metrics])),
        "auprc_oof": metric_value(y_true, scores, "auprc"),
        "auprc_ci_low": auprc_ci[0],
        "auprc_ci_high": auprc_ci[1],
        "auprc_fold_mean": float(np.nanmean([item["auprc"] for item in fold_metrics])),
        "auprc_fold_std": float(np.nanstd([item["auprc"] for item in fold_metrics])),
        "brier_oof": metric_value(y_true, scores, "brier"),
        "brier_ci_low": brier_ci[0],
        "brier_ci_high": brier_ci[1],
        "brier_fold_mean": float(np.nanmean([item["brier"] for item in fold_metrics])),
        "brier_fold_std": float(np.nanstd([item["brier"] for item in fold_metrics])),
    }


def evaluate_feature_set(
    df: pd.DataFrame,
    y: np.ndarray,
    columns: list[str],
    name: str,
    splits: list[tuple[np.ndarray, np.ndarray]],
    args: argparse.Namespace,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    x = numeric_feature_frame(df, columns)
    predictions: list[dict[str, object]] = []
    fold_metrics = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        model = make_model(args.model)
        model.fit(x.iloc[train_idx], y[train_idx])
        scores = model_scores(model, x.iloc[test_idx])
        y_test = y[test_idx]
        fold_metrics.append(
            {
                "auroc": metric_value(y_test, scores, "auroc"),
                "auprc": metric_value(y_test, scores, "auprc"),
                "brier": metric_value(y_test, scores, "brier"),
            }
        )
        for row_index, score in zip(test_idx, scores):
            predictions.append(
                {
                    "feature_set": name,
                    "fold": fold,
                    "variant_id": df.iloc[row_index].get("variant_id", ""),
                    "gene_id": df.iloc[row_index].get("gene_id", ""),
                    "strategy": df.iloc[row_index].get("strategy", ""),
                    "label": int(y[row_index]),
                    "score": float(score),
                }
            )
    pred_df = pd.DataFrame(predictions)
    return (
        {
            "feature_set": name,
            "feature_count": len(columns),
            "row_count": len(df),
            **summarize_predictions(pred_df, fold_metrics, args),
        },
        predictions,
    )


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.dataset_tsv, sep="\t", low_memory=False)
    if args.label_column not in df.columns:
        raise ValueError(f"Missing label column: {args.label_column}")
    labels = df[args.label_column].map(parse_label)
    df = df[labels.notna()].copy()
    y = labels[labels.notna()].astype(int).to_numpy()
    if len(df) == 0:
        raise ValueError("No rows with usable binary labels")
    if len(set(y)) != 2:
        raise ValueError("Evaluation requires both positive and negative labels")

    gaph_features = infer_gaph_features(df, args.gaph_prefix)
    if not gaph_features:
        raise ValueError(f"No numeric GAPH features found with prefix {args.gaph_prefix!r}")
    baseline_features = parse_feature_list(args.baseline_features)
    if not baseline_features:
        baseline_features = infer_baseline_features(df, args.label_column, gaph_features)
    if not baseline_features:
        raise ValueError("No baseline features supplied or inferred")

    shuffled = df.copy()
    rng = np.random.default_rng(args.random_seed)
    for col in gaph_features:
        shuffled[col] = rng.permutation(shuffled[col].to_numpy())

    feature_sets = [
        ("baseline", df, baseline_features),
        ("gaph", df, gaph_features),
        ("baseline_plus_gaph", df, baseline_features + gaph_features),
        ("baseline_plus_shuffled_gaph", shuffled, baseline_features + gaph_features),
    ]

    args.outdir.mkdir(parents=True, exist_ok=True)
    splits = split_indices(df, y, args.group_column, args.folds, args.random_seed)
    split_rows = split_assignment_rows(df, splits)
    write_tsv(args.outdir / "split_assignments.tsv", split_rows, list(split_rows[0].keys()))
    metric_rows = []
    prediction_rows = []
    for name, frame, columns in feature_sets:
        metrics, predictions = evaluate_feature_set(frame, y, columns, name, splits, args)
        metric_rows.append(metrics)
        prediction_rows.extend(predictions)

    write_tsv(args.outdir / "metrics.tsv", metric_rows, list(metric_rows[0].keys()))
    write_tsv(args.outdir / "predictions.tsv", prediction_rows, list(prediction_rows[0].keys()))
    summary = {
        "dataset_tsv": str(args.dataset_tsv),
        "row_count": int(len(df)),
        "positive_count": int(y.sum()),
        "negative_count": int(len(y) - y.sum()),
        "baseline_features": baseline_features,
        "gaph_features": gaph_features,
        "model": args.model,
        "group_column": args.group_column,
        "folds": args.folds,
        "bootstrap_iterations": args.bootstrap_iterations,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
