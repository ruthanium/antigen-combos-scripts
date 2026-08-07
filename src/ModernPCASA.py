"""Patient-aware, interpretable PCASA-style antigen gate optimization.

Discovery uses tumor and normal epithelial cells from GSE178341. Entire
patients, rather than cells, are assigned to train or test. Top gates are then
validated against explicitly annotated cancer cells and normal epithelial
cells in the independent ICBI CRC atlas.
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
from pathlib import Path
import time
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

from ValidateColonSingleCell import (
    get_icbi_gene,
    get_icbi_metadata,
    read_hacohen_target_counts,
)


METADATA = {"dataset", "batch", "type", "tissue.cancer"}
GENE_ALIASES = {"PVRL1": "NECTIN1", "GPR133": "ADGRD1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("data/test-normalized-matrix.txt"))
    parser.add_argument("--data-dir", type=Path, default=Path("data/single_cell"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/modern_pcasa")
    )
    parser.add_argument("--top-genes", type=int, default=30)
    parser.add_argument("--top-per-logic", type=int, default=25)
    parser.add_argument("--external-per-logic", type=int, default=10)
    parser.add_argument("--test-size", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--external-only",
        action="store_true",
        help="Resume ICBI validation from selected_for_external_validation.csv",
    )
    return parser.parse_args()


def candidate_genes(matrix_path: Path) -> list[str]:
    columns = pd.read_csv(matrix_path, nrows=0).columns
    return [
        GENE_ALIASES.get(column, column)
        for column in columns
        if column not in METADATA
    ]


def patient_rates(
    predictions: np.ndarray, labels: np.ndarray, patients: np.ndarray, positive: bool
) -> np.ndarray:
    frame = pd.DataFrame(
        {
            "prediction": predictions[labels == positive],
            "patient": patients[labels == positive],
        }
    )
    if frame.empty:
        return np.array([], dtype=float)
    return frame.groupby("patient", observed=True)["prediction"].mean().to_numpy()


def gate_metrics(
    predictions: np.ndarray,
    labels: np.ndarray,
    patients: np.ndarray,
    prefix: str,
) -> dict[str, float]:
    predictions = predictions.astype(bool, copy=False)
    labels = labels.astype(bool, copy=False)
    tp = int(np.sum(predictions & labels))
    fp = int(np.sum(predictions & ~labels))
    fn = int(np.sum(~predictions & labels))
    tn = int(np.sum(~predictions & ~labels))
    tumor_coverage = tp / (tp + fn) if tp + fn else np.nan
    normal_coverage = fp / (fp + tn) if fp + tn else np.nan
    precision = tp / (tp + fp) if tp + fp else np.nan
    f1 = (
        2 * precision * tumor_coverage / (precision + tumor_coverage)
        if np.isfinite(precision) and precision + tumor_coverage
        else np.nan
    )
    tumor_patient = patient_rates(predictions, labels, patients, True)
    normal_patient = patient_rates(predictions, labels, patients, False)
    patient_prevalence = (
        float(np.mean(tumor_patient >= 0.01)) if tumor_patient.size else np.nan
    )
    # Safety dominates progressively as normal coverage rises. Patient
    # prevalence prevents a gate driven by a small number of exceptional cases.
    score = (
        tumor_coverage * (1 - normal_coverage) ** 5 * patient_prevalence
        if np.isfinite(tumor_coverage)
        and np.isfinite(normal_coverage)
        and np.isfinite(patient_prevalence)
        else np.nan
    )
    return {
        f"{prefix}_tumor_coverage": tumor_coverage,
        f"{prefix}_normal_coverage": normal_coverage,
        f"{prefix}_precision": precision,
        f"{prefix}_f1": f1,
        f"{prefix}_patient_prevalence": patient_prevalence,
        f"{prefix}_median_patient_tumor_coverage": (
            float(np.median(tumor_patient)) if tumor_patient.size else np.nan
        ),
        f"{prefix}_max_patient_normal_coverage": (
            float(np.max(normal_patient)) if normal_patient.size else np.nan
        ),
        f"{prefix}_score": score,
    }


def prediction_for_gate(
    matrix: np.ndarray, indices: Sequence[int], logic: str
) -> np.ndarray:
    selected = matrix[:, indices]
    if logic == "SINGLE":
        return selected[:, 0]
    if logic == "AND":
        return np.all(selected, axis=1)
    if logic == "OR":
        return np.any(selected, axis=1)
    if logic == "NOT":
        return selected[:, 0] & ~selected[:, 1]
    raise ValueError(f"Unknown logic: {logic}")


def evaluate_candidates(
    expression: np.ndarray,
    labels: np.ndarray,
    patients: np.ndarray,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    genes: Sequence[str],
    gates: Iterable[Sequence[int]],
    logic: str,
) -> pd.DataFrame:
    rows = []
    for indices in gates:
        indices = tuple(indices)
        prediction = prediction_for_gate(expression, indices, logic)
        row: dict[str, object] = {
            "logic": logic,
            "gate": (
                f"{genes[indices[0]]} NOT {genes[indices[1]]}"
                if logic == "NOT"
                else f" {logic} ".join(genes[index] for index in indices)
            ),
            "genes": "|".join(genes[index] for index in indices),
        }
        row.update(
            gate_metrics(
                prediction[train_mask],
                labels[train_mask],
                patients[train_mask],
                "train",
            )
        )
        row.update(
            gate_metrics(
                prediction[test_mask],
                labels[test_mask],
                patients[test_mask],
                "test",
            )
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["train_score", "train_tumor_coverage"],
        ascending=False,
        kind="stable",
    )


def load_hacohen(
    data_dir: Path, genes: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    stem = "GSE178341_crc10x_full_c295v4_submit"
    expression, barcodes = read_hacohen_target_counts(
        data_dir / f"{stem}.h5", genes, 10_000
    )
    clusters = pd.read_csv(data_dir / f"{stem}_cluster.csv.gz")
    clinical = pd.read_csv(data_dir / f"{stem}_metatables.csv.gz")
    if not np.array_equal(barcodes, clusters["sampleID"].to_numpy()):
        raise ValueError("Hacohen expression and clusters are not aligned")
    epithelial = clusters["clTopLevel"].eq("Epi").to_numpy()
    specimen = clinical["SPECIMEN_TYPE"].to_numpy()
    keep = epithelial & np.isin(specimen, ["T", "N"])
    detected = expression.to_numpy()[keep] > 0
    labels = specimen[keep] == "T"
    patients = clinical.loc[keep, "PID"].astype(str).to_numpy()
    return detected, labels, patients


def patient_split(
    labels: np.ndarray, patients: np.ndarray, test_size: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
    train_index, test_index = next(
        splitter.split(np.zeros(len(labels)), labels, groups=patients)
    )
    train = np.zeros(len(labels), dtype=bool)
    test = np.zeros(len(labels), dtype=bool)
    train[train_index] = True
    test[test_index] = True
    return train, test


def select_external_gates(
    rankings: Sequence[pd.DataFrame], per_logic: int
) -> pd.DataFrame:
    selected = []
    for ranking in rankings:
        # Selection is based on training metrics. Include both the overall
        # utility leaders and explicit safety-frontier candidates.
        stable = ranking.loc[ranking["test_tumor_coverage"].gt(0)]
        selected.append(stable.head(per_logic))
        for ceiling in (0.01, 0.001):
            safe = stable.loc[stable["train_normal_coverage"].le(ceiling)].sort_values(
                ["train_tumor_coverage", "train_patient_prevalence"],
                ascending=False,
                kind="stable",
            )
            selected.append(safe.head(per_logic))
    return (
        pd.concat(selected, ignore_index=True)
        .drop_duplicates(["logic", "gate"])
        .reset_index(drop=True)
    )


def load_icbi_cached(
    genes: Sequence[str], cache_dir: Path
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = cache_dir / "metadata.pkl"
    if metadata_path.exists():
        metadata = pd.read_pickle(metadata_path)
    else:
        metadata = None
    if metadata is None or "dataset" not in metadata.columns:
        metadata = get_icbi_metadata()
        metadata.to_pickle(metadata_path)

    columns = []
    available_genes = []
    for gene in genes:
        path = cache_dir / f"{gene}.npy"
        if path.exists():
            values = np.load(path)
        else:
            print(f"ICBI external validation: fetching {gene}", flush=True)
            values = None
            for attempt in range(1, 5):
                try:
                    values = get_icbi_gene(gene)
                    break
                except Exception as error:
                    if attempt == 4:
                        print(
                            f"Skipping {gene}: ICBI transfer failed after 4 attempts",
                            flush=True,
                        )
                        break
                    print(
                        f"Retrying {gene} after transfer failure ({attempt}/4)",
                        flush=True,
                    )
                    time.sleep(2**attempt)
            if values is None:
                continue
            np.save(path, values)
        available_genes.append(gene)
        columns.append(values > 0)
    return metadata, np.column_stack(columns), available_genes


def external_validate(
    selected: pd.DataFrame, output_dir: Path
) -> pd.DataFrame:
    genes = sorted(
        {gene for value in selected["genes"] for gene in str(value).split("|")}
    )
    metadata, detected, available_genes = load_icbi_cached(
        genes, output_dir / "icbi_cache"
    )
    positive = metadata["sample_type"].eq("primary tumor") & metadata[
        "atlas_cell_type_coarse"
    ].eq("Cancer cell")
    negative = metadata["sample_type"].isin(
        ["adjacent normal", "healthy normal"]
    ) & metadata["atlas_cell_type_coarse"].eq("Epithelial cell")
    discovery_overlap = metadata["dataset"].astype(str).str.contains(
        "Pelka|GSE178341", case=False, regex=True
    )
    keep = ((positive | negative) & ~discovery_overlap).to_numpy()
    labels = positive.to_numpy()[keep]
    patients = metadata.loc[keep, "patient_id"].astype(str).to_numpy()
    detected = detected[keep]
    gene_index = {gene: index for index, gene in enumerate(available_genes)}

    rows = []
    for _, gate in selected.iterrows():
        gate_genes = str(gate["genes"]).split("|")
        if not set(gate_genes).issubset(gene_index):
            continue
        indices = [gene_index[gene] for gene in gate_genes]
        prediction = prediction_for_gate(detected, indices, str(gate["logic"]))
        row = gate.to_dict()
        row.update(gate_metrics(prediction, labels, patients, "external"))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["external_score", "external_tumor_coverage"],
        ascending=False,
        kind="stable",
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.external_only:
        selected = pd.read_csv(
            args.output_dir / "selected_for_external_validation.csv"
        )
        external = external_validate(selected, args.output_dir)
        external.to_csv(args.output_dir / "external_validation.csv", index=False)
        print(f"External results written to {args.output_dir}", flush=True)
        return

    genes = candidate_genes(args.matrix)
    detected, labels, patients = load_hacohen(args.data_dir, genes)
    train, test = patient_split(labels, patients, args.test_size, args.seed)

    split_summary = pd.DataFrame(
        {
            "split": ["train", "test"],
            "cells": [int(train.sum()), int(test.sum())],
            "tumor_cells": [int(np.sum(labels & train)), int(np.sum(labels & test))],
            "normal_cells": [int(np.sum(~labels & train)), int(np.sum(~labels & test))],
            "patients": [
                len(np.unique(patients[train])),
                len(np.unique(patients[test])),
            ],
        }
    )
    split_summary.to_csv(args.output_dir / "discovery_split.csv", index=False)

    print(f"Evaluating {len(genes)} single genes", flush=True)
    singles = evaluate_candidates(
        detected,
        labels,
        patients,
        train,
        test,
        genes,
        ((index,) for index in range(len(genes))),
        "SINGLE",
    )
    singles.to_csv(args.output_dir / "single_rankings.csv", index=False)
    top_gene_names = (
        singles.loc[singles["train_tumor_coverage"].ge(0.01), "genes"]
        .head(args.top_genes)
        .tolist()
    )
    top_indices = [genes.index(gene) for gene in top_gene_names]

    print(f"Evaluating gates among {len(top_indices)} genes", flush=True)
    and_pairs = evaluate_candidates(
        detected,
        labels,
        patients,
        train,
        test,
        genes,
        combinations(top_indices, 2),
        "AND",
    )
    or_pairs = evaluate_candidates(
        detected,
        labels,
        patients,
        train,
        test,
        genes,
        combinations(top_indices, 2),
        "OR",
    )
    not_pairs = evaluate_candidates(
        detected,
        labels,
        patients,
        train,
        test,
        genes,
        permutations(top_indices, 2),
        "NOT",
    )
    and_triples = evaluate_candidates(
        detected,
        labels,
        patients,
        train,
        test,
        genes,
        combinations(top_indices, 3),
        "AND",
    )
    rankings = [singles, and_pairs, or_pairs, not_pairs, and_triples]
    names = ["single", "and_pairs", "or_pairs", "not_pairs", "and_triples"]
    for name, ranking in zip(names, rankings):
        ranking.to_csv(args.output_dir / f"all_{name}.csv", index=False)
        ranking.head(args.top_per_logic).to_csv(
            args.output_dir / f"top_{name}.csv", index=False
        )

    selected = select_external_gates(rankings, args.external_per_logic)
    selected.to_csv(args.output_dir / "selected_for_external_validation.csv", index=False)
    external = external_validate(selected, args.output_dir)
    external.to_csv(args.output_dir / "external_validation.csv", index=False)
    print(f"Results written to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
