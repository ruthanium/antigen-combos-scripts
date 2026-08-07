"""Shared functions for the Python port of the antigen-combination pipeline."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier


METADATA_COLUMNS = {"dataset", "batch", "type", "tissue.cancer"}


def read_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = METADATA_COLUMNS
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Expression matrix is missing columns: {sorted(missing)}")
    return frame.drop(columns=["batch"])


def gene_columns(frame: pd.DataFrame) -> list[str]:
    return [column for column in frame.columns if column not in METADATA_COLUMNS]


def cancer_types(frame: pd.DataFrame) -> list[str]:
    return frame.loc[frame["type"].eq("cancer"), "tissue.cancer"].drop_duplicates().tolist()


def sketch_path(results_dir: Path, cancer: str) -> Path:
    name = cancer.lower().replace(" ", "-") + "-sketch.txt"
    return results_dir / "sketches" / name


def cancer_normal_subset(
    frame: pd.DataFrame, cancer: str, sample_ids: set[str] | None = None
) -> pd.DataFrame:
    selected = frame["type"].eq("tissue") | frame["tissue.cancer"].eq(cancer)
    if sample_ids is not None:
        selected &= frame["dataset"].isin(sample_ids)
    return frame.loc[selected]


def read_sketch(results_dir: Path, cancer: str) -> set[str]:
    path = sketch_path(results_dir, cancer)
    if not path.exists():
        raise FileNotFoundError(f"Missing sketch for {cancer}: {path}")
    return set(path.read_text(encoding="utf-8").splitlines())


def _r_fixed_db(values: np.ndarray, target_mask: np.ndarray) -> tuple[float, float, str]:
    """Reproduce index.DB.fixed(..., centrotypes='centroids', p=1, q=2).

    The original function intentionally differs from the conventional
    Davies–Bouldin implementation in its within-cluster scatter formula.
    """
    target = values[target_mask]
    other = values[~target_mask]
    if len(target) == 0 or len(other) == 0:
        return np.nan, np.nan, ""

    target_center = target.mean(axis=0)
    other_center = other.mean(axis=0)

    def scatter(group: np.ndarray, center: np.ndarray) -> float:
        if len(group) <= 1:
            return 0.0
        squared_distance = np.square(group - center).sum(axis=1)
        return float(np.sqrt(np.mean(np.square(squared_distance))))

    separation = float(np.abs(target_center - other_center).sum())
    if separation == 0:
        db = np.inf
    else:
        db = (scatter(target, target_center) + scatter(other, other_center)) / separation

    gate_type = "".join(
        "H" if target_value >= other_value else "L"
        for target_value, other_value in zip(target_center, other_center)
    )
    return db, separation, gate_type


def clustering_scores(
    frame: pd.DataFrame,
    cancer: str,
    combos: Iterable[Sequence[str]],
) -> pd.DataFrame:
    target_mask = frame["tissue.cancer"].eq(cancer).to_numpy()
    rows: list[dict[str, object]] = []
    for combo in combos:
        genes = tuple(combo)
        values = frame.loc[:, list(genes)].to_numpy(dtype=float)
        db, distance, gate_type = _r_fixed_db(values, target_mask)
        rows.append(
            {
                "can": cancer,
                "combo": ":".join(genes),
                "db": db,
                "dist.man": distance,
                "ttype": gate_type,
            }
        )
    return pd.DataFrame(rows)


def all_combinations(genes: Sequence[str], size: int, limit: int | None = None):
    iterator = combinations(genes, size)
    if limit is None:
        return iterator
    from itertools import islice

    return islice(iterator, limit)


def add_antigen_labels(
    scores: pd.DataFrame, clinical: set[str], size: int
) -> pd.DataFrame:
    result = scores.copy()
    split = result["combo"].str.split(":", expand=True)
    labels = []
    for index in range(size):
        gene_column = f"gene{index + 1}"
        label_column = f"g{index + 1}lab"
        result[gene_column] = split[index]
        result[label_column] = np.where(result[gene_column].isin(clinical), "C", "N")
        labels.append(result[label_column])
    result["lab"] = labels[0]
    for label in labels[1:]:
        result["lab"] += ":" + label
    return result


def _fit_stump_predictions(
    train_values: np.ndarray,
    train_target: np.ndarray,
    test_values: np.ndarray,
) -> np.ndarray:
    """Fit the closest sklearn equivalent of the rpart depth-one classifier."""
    if len(np.unique(train_target)) < 2:
        return np.full(len(test_values), bool(train_target[0]), dtype=bool)
    model = DecisionTreeClassifier(
        criterion="gini",
        max_depth=1,
        min_samples_split=20,
        min_samples_leaf=7,
        random_state=0,
    )
    model.fit(train_values.reshape(-1, 1), train_target)
    return model.predict(test_values.reshape(-1, 1)).astype(bool)


def evaluate_combo(
    frame: pd.DataFrame,
    cancer: str,
    genes: Sequence[str],
    sketch_ids: set[str],
) -> dict[str, float]:
    relevant = cancer_normal_subset(frame, cancer)
    is_train = relevant["dataset"].isin(sketch_ids).to_numpy()
    train = relevant.loc[is_train]
    test = relevant.loc[~is_train]
    if train.empty or test.empty:
        return {"f1": np.nan, "prec": np.nan, "rec": np.nan}

    train_target = train["tissue.cancer"].eq(cancer).to_numpy()
    test_target = test["tissue.cancer"].eq(cancer).to_numpy()
    gene_predictions = []
    for gene in genes:
        gene_predictions.append(
            _fit_stump_predictions(
                train[gene].to_numpy(dtype=float),
                train_target,
                test[gene].to_numpy(dtype=float),
            )
        )

    # pmin over {-1, 1} in the R code is an AND gate.
    predicted_target = np.logical_and.reduce(gene_predictions)
    tp = int(np.sum(test_target & predicted_target))
    fp = int(np.sum(~test_target & predicted_target))
    positives = int(np.sum(test_target))
    precision = tp / (tp + fp) if tp + fp else np.nan
    recall = tp / positives if positives else np.nan
    f1 = (
        2 * precision * recall / (precision + recall)
        if np.isfinite(precision) and precision + recall
        else np.nan
    )
    return {"f1": f1, "prec": precision, "rec": recall}


def clean_and_scale(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.loc[frame["db"].notna()].copy()

    def scale(series: pd.Series, reverse: bool = False) -> pd.Series:
        logged = np.log(series.astype(float))
        span = logged.max() - logged.min()
        if not np.isfinite(span) or span == 0:
            scaled = pd.Series(1.0, index=series.index)
        else:
            scaled = (logged - logged.min()) / span
        return 1 - scaled if reverse else scaled

    result["dbscaled"] = scale(result["db"], reverse=True)
    result["distscaled"] = scale(result["dist.man"])
    result["score"] = result[["dbscaled", "distscaled"]].min(axis=1)
    return result.sort_values("score", ascending=False, kind="stable")


def diverse_top_k(
    frame: pd.DataFrame, size: int, k: int = 10, max_gene_uses: int = 2
) -> pd.DataFrame:
    chosen = []
    counts: dict[str, int] = {}
    for _, row in frame.iterrows():
        genes = row["combo"].split(":")
        if any(counts.get(gene, 0) >= max_gene_uses for gene in genes):
            continue
        chosen.append(row)
        for gene in genes:
            counts[gene] = counts.get(gene, 0) + 1
        if len(chosen) >= k:
            break
    if not chosen:
        return frame.iloc[0:0].copy()
    return pd.DataFrame(chosen).reset_index(drop=True)
