"""Find two single-cell RNA partners for an anchored protein biomarker.

The anchor can be supplied as measured protein-positive cell barcodes (for
example, CITE-seq/flow-derived calls). When protein calls are unavailable,
detected anchor RNA is used and clearly reported as an RNA surrogate.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

from ModernPCASA import candidate_genes, load_icbi_cached
from ValidateColonSingleCell import read_hacohen_target_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--biomarker", required=True)
    parser.add_argument("--extra-genes", nargs="*", default=[])
    parser.add_argument("--candidate-file", type=Path)
    parser.add_argument("--protein-calls", type=Path)
    parser.add_argument("--protein-threshold", type=float, default=0.0)
    parser.add_argument("--top-partners", type=int, default=30)
    parser.add_argument("--external-top", type=int, default=10)
    parser.add_argument("--max-normal", type=float, default=0.01)
    parser.add_argument("--test-size", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--data-dir", type=Path, default=Path("data/single_cell"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/biomarker_partners")
    )
    return parser.parse_args()


def phi(left: np.ndarray, right: np.ndarray) -> float:
    left = left.astype(float, copy=False)
    right = right.astype(float, copy=False)
    if left.std() == 0 or right.std() == 0:
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def patient_summary(
    positive: np.ndarray, patients: np.ndarray
) -> tuple[float, float, float]:
    rates = (
        pd.DataFrame({"positive": positive, "patient": patients})
        .groupby("patient", observed=True)["positive"]
        .mean()
    )
    return float((rates >= 0.01).mean()), float(rates.median()), float(rates.max())


def gate_metrics(
    positive: np.ndarray,
    labels: np.ndarray,
    patients: np.ndarray,
) -> dict[str, float]:
    row: dict[str, float] = {}
    for value, name in ((True, "tumor"), (False, "normal")):
        selected = labels == value
        calls = positive[selected]
        prevalence, median, maximum = patient_summary(calls, patients[selected])
        row[f"{name}_coverage"] = float(calls.mean())
        row[f"{name}_patient_prevalence"] = prevalence
        row[f"{name}_median_patient_coverage"] = median
        row[f"{name}_max_patient_coverage"] = maximum
    row["score"] = (
        row["tumor_coverage"]
        * (1 - row["normal_coverage"]) ** 8
        * row["tumor_patient_prevalence"]
    )
    return row


def prefixed_metrics(
    positive: np.ndarray,
    labels: np.ndarray,
    patients: np.ndarray,
    mask: np.ndarray,
    prefix: str,
) -> dict[str, float]:
    return {
        f"{prefix}_{key}": value
        for key, value in gate_metrics(
            positive[mask], labels[mask], patients[mask]
        ).items()
    }


def bootstrap_patient_ci(
    positive: np.ndarray,
    labels: np.ndarray,
    patients: np.ndarray,
    mask: np.ndarray,
    repeats: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    result: dict[str, float] = {}
    for label, name in ((True, "tumor"), (False, "normal")):
        selected = mask & (labels == label)
        rates = (
            pd.DataFrame(
                {"positive": positive[selected], "patient": patients[selected]}
            )
            .groupby("patient", observed=True)["positive"]
            .mean()
            .to_numpy()
        )
        samples = np.array(
            [rng.choice(rates, len(rates), replace=True).mean() for _ in range(repeats)]
        )
        result[f"test_{name}_patient_mean_ci_low"] = float(
            np.quantile(samples, 0.025)
        )
        result[f"test_{name}_patient_mean_ci_high"] = float(
            np.quantile(samples, 0.975)
        )
    return result


def read_protein_anchor(
    path: Path,
    barcodes: np.ndarray,
    threshold: float,
) -> np.ndarray:
    calls = pd.read_csv(path, sep=None, engine="python")
    barcode_column = next(
        (column for column in ("cell_id", "barcode", "cellID") if column in calls),
        None,
    )
    if barcode_column is None:
        raise ValueError("Protein calls need a cell_id, barcode, or cellID column")
    value_columns = [column for column in calls if column != barcode_column]
    if value_columns:
        positive_ids = calls.loc[
            pd.to_numeric(calls[value_columns[0]], errors="coerce") > threshold,
            barcode_column,
        ]
    else:
        positive_ids = calls[barcode_column]
    return np.isin(barcodes, positive_ids.astype(str).to_numpy())


def cluster_decode(
    frame: pd.DataFrame,
    calls: np.ndarray,
    gate: str,
) -> pd.DataFrame:
    selected = frame.assign(gate_positive=calls)
    grouped = (
        selected.groupby(
            ["SPECIMEN_TYPE", "cl295v11SubShort", "cl295v11SubFull"],
            observed=True,
        )
        .agg(
            cells=("gate_positive", "size"),
            positive_cells=("gate_positive", "sum"),
            coverage=("gate_positive", "mean"),
            patients=("PID", "nunique"),
        )
        .reset_index()
    )
    grouped.insert(0, "gate", gate)
    return grouped


def marker_expression_summary(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    group_columns: list[str],
    source: str,
    expression_units: str,
) -> pd.DataFrame:
    joined = pd.concat(
        [metadata.reset_index(drop=True), expression[genes].reset_index(drop=True)],
        axis=1,
    )
    rows = []
    for group_values, group in joined.groupby(
        group_columns, observed=True, dropna=False
    ):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        group_info = dict(zip(group_columns, group_values))
        for gene in genes:
            values = group[gene].to_numpy(dtype=float)
            positive = values > 0
            rows.append(
                {
                    "source": source,
                    "expression_units": expression_units,
                    **group_info,
                    "marker": gene,
                    "cells": len(values),
                    "positive_cells": int(positive.sum()),
                    "pct_expressed": 100 * float(positive.mean()),
                    "mean_expression_all_cells": float(values.mean()),
                    "median_expression_all_cells": float(np.median(values)),
                    "mean_expression_positive_cells": (
                        float(values[positive].mean()) if positive.any() else 0.0
                    ),
                    "median_expression_positive_cells": (
                        float(np.median(values[positive])) if positive.any() else 0.0
                    ),
                    "q25_expression_positive_cells": (
                        float(np.quantile(values[positive], 0.25))
                        if positive.any()
                        else 0.0
                    ),
                    "q75_expression_positive_cells": (
                        float(np.quantile(values[positive], 0.75))
                        if positive.any()
                        else 0.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def external_validate(
    triples: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    selected_genes = sorted(
        {
            gene
            for _, row in triples.iterrows()
            for gene in (row["biomarker"], row["partner_1"], row["partner_2"])
        }
    )
    metadata, detected, available = load_icbi_cached(
        selected_genes, output_dir / "icbi_cache"
    )
    positive = metadata["sample_type"].eq("primary tumor") & metadata[
        "atlas_cell_type_coarse"
    ].eq("Cancer cell")
    negative = metadata["sample_type"].isin(
        ["adjacent normal", "healthy normal"]
    ) & metadata["atlas_cell_type_coarse"].eq("Epithelial cell")
    overlap = metadata["dataset"].astype(str).str.contains(
        "Pelka|GSE178341", case=False, regex=True
    )
    keep = ((positive | negative) & ~overlap).to_numpy()
    labels = positive.to_numpy()[keep]
    patients = metadata.loc[keep, "patient_id"].astype(str).to_numpy()
    detected = detected[keep]
    index = {gene: position for position, gene in enumerate(available)}
    raw_expression = pd.DataFrame(
        {
            gene: np.load(output_dir / "icbi_cache" / f"{gene}.npy")[keep]
            for gene in available
        }
    )
    external_metadata = pd.DataFrame(
        {
            "compartment": np.where(labels, "cancer_cell", "normal_epithelial"),
            "patient_id": patients,
        }
    )
    marker_expression_summary(
        raw_expression,
        external_metadata,
        available,
        ["compartment"],
        "ICBI_CRC_atlas_nonoverlap",
        "ICBI_atlas_expression_value",
    ).to_csv(output_dir / "external_marker_expression.csv", index=False)
    rows = []
    for _, gate in triples.iterrows():
        genes = [gate["biomarker"], gate["partner_1"], gate["partner_2"]]
        if not set(genes).issubset(index):
            continue
        calls = np.all(detected[:, [index[gene] for gene in genes]], axis=1)
        row = gate.to_dict()
        row.update(
            {
                f"external_{key}": value
                for key, value in gate_metrics(calls, labels, patients).items()
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["external_score", "external_tumor_coverage"],
        ascending=False,
        kind="stable",
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.candidate_file:
        candidates = pd.read_csv(args.candidate_file, header=None)[0].astype(str).tolist()
    else:
        candidates = candidate_genes(Path("data/test-normalized-matrix.txt"))
    genes = list(dict.fromkeys([args.biomarker, *candidates, *args.extra_genes]))

    stem = "GSE178341_crc10x_full_c295v4_submit"
    expression, barcodes = read_hacohen_target_counts(
        args.data_dir / f"{stem}.h5", genes, 10_000
    )
    clusters = pd.read_csv(args.data_dir / f"{stem}_cluster.csv.gz")
    clinical = pd.read_csv(args.data_dir / f"{stem}_metatables.csv.gz")
    if not np.array_equal(barcodes, clusters["sampleID"].to_numpy()):
        raise ValueError("Expression and cluster annotations are not aligned")
    frame = pd.concat(
        [clusters, clinical.drop(columns=["cellID"]), expression.gt(0)],
        axis=1,
    )
    keep = frame["clTopLevel"].eq("Epi") & frame["SPECIMEN_TYPE"].isin(["T", "N"])
    frame = frame.loc[keep].reset_index(drop=True)
    detected = frame[genes].to_numpy(dtype=bool)
    labels = frame["SPECIMEN_TYPE"].eq("T").to_numpy()
    patients = frame["PID"].astype(str).to_numpy()
    splitter = GroupShuffleSplit(
        n_splits=1, test_size=args.test_size, random_state=args.seed
    )
    train_index, test_index = next(
        splitter.split(np.zeros(len(labels)), labels, groups=patients)
    )
    train = np.zeros(len(labels), dtype=bool)
    test = np.zeros(len(labels), dtype=bool)
    train[train_index] = True
    test[test_index] = True

    if args.protein_calls:
        full_anchor = read_protein_anchor(
            args.protein_calls, barcodes.astype(str), args.protein_threshold
        )
        anchor = full_anchor[keep.to_numpy()]
        anchor_source = "measured_protein"
    else:
        anchor = detected[:, genes.index(args.biomarker)]
        anchor_source = "detected_RNA_surrogate"

    partner_rows = []
    for index, gene in enumerate(genes):
        if gene == args.biomarker:
            continue
        partner = detected[:, index]
        row: dict[str, object] = {
            "biomarker": args.biomarker,
            "anchor_source": anchor_source,
            "partner": gene,
            "tumor_phi": phi(anchor[labels], partner[labels]),
            "normal_phi": phi(anchor[~labels], partner[~labels]),
            "tumor_partner_given_anchor": float(partner[labels & anchor].mean()),
            "normal_partner_given_anchor": float(partner[~labels & anchor].mean()),
        }
        calls = anchor & partner
        row.update(prefixed_metrics(calls, labels, patients, train, "train"))
        row.update(prefixed_metrics(calls, labels, patients, test, "test"))
        partner_rows.append(row)
    partner_table = pd.DataFrame(partner_rows).sort_values(
        ["train_score", "train_tumor_coverage"], ascending=False, kind="stable"
    )
    partner_table.to_csv(args.output_dir / "partner_ranking.csv", index=False)

    eligible = partner_table[
        partner_table["train_normal_coverage"] <= args.max_normal
    ].head(args.top_partners)
    eligible_genes = eligible["partner"].tolist()
    index_by_gene = {gene: genes.index(gene) for gene in genes}
    triple_rows = []
    cluster_tables = []
    for left, right in combinations(eligible_genes, 2):
        calls = (
            anchor
            & detected[:, index_by_gene[left]]
            & detected[:, index_by_gene[right]]
        )
        gate = f"{args.biomarker} AND {left} AND {right}"
        row: dict[str, object] = {
            "biomarker": args.biomarker,
            "anchor_source": anchor_source,
            "partner_1": left,
            "partner_2": right,
            "gate": gate,
        }
        row.update(prefixed_metrics(calls, labels, patients, train, "train"))
        row.update(prefixed_metrics(calls, labels, patients, test, "test"))
        row.update(
            bootstrap_patient_ci(
                calls, labels, patients, test, args.bootstrap, args.seed
            )
        )
        triple_rows.append(row)
        cluster_tables.append(
            cluster_decode(frame.loc[test].reset_index(drop=True), calls[test], gate)
        )
    triples = pd.DataFrame(triple_rows).sort_values(
        ["train_score", "train_tumor_coverage"], ascending=False, kind="stable"
    )
    triples.to_csv(args.output_dir / "two_partner_gates.csv", index=False)
    if cluster_tables:
        pd.concat(cluster_tables, ignore_index=True).to_csv(
            args.output_dir / "two_partner_cluster_decode.csv", index=False
        )
    external_candidates = (
        triples.loc[triples["test_tumor_coverage"].gt(0)]
        .sort_values(
            ["test_score", "test_tumor_coverage"],
            ascending=False,
            kind="stable",
        )
        .head(args.external_top)
    )
    external_candidates.to_csv(
        args.output_dir / "selected_for_external_validation.csv", index=False
    )
    selected_markers = sorted(
        {
            gene
            for _, row in external_candidates.iterrows()
            for gene in (row["biomarker"], row["partner_1"], row["partner_2"])
        }
    )
    split_metadata = frame[
        ["SPECIMEN_TYPE", "PID", "cl295v11SubShort", "cl295v11SubFull"]
    ].copy()
    split_metadata["split"] = np.where(train, "discovery", "held_out")
    split_metadata["compartment"] = np.where(
        labels, "tumor_epithelial", "normal_epithelial"
    )
    marker_expression_summary(
        expression.loc[keep].reset_index(drop=True),
        split_metadata,
        selected_markers,
        ["split", "compartment"],
        "GSE178341",
        "raw_UMI_counts",
    ).to_csv(args.output_dir / "marker_expression_by_split.csv", index=False)
    heldout = test
    marker_expression_summary(
        expression.loc[keep].reset_index(drop=True).loc[heldout].reset_index(drop=True),
        split_metadata.loc[heldout].reset_index(drop=True),
        selected_markers,
        ["compartment", "cl295v11SubShort", "cl295v11SubFull"],
        "GSE178341_held_out",
        "raw_UMI_counts",
    ).to_csv(
        args.output_dir / "heldout_marker_expression_by_cluster.csv", index=False
    )
    external = external_validate(external_candidates, args.output_dir)
    external.to_csv(args.output_dir / "external_validation.csv", index=False)
    pd.DataFrame(
        {
            "biomarker": [args.biomarker],
            "anchor_source": [anchor_source],
            "epithelial_cells": [len(frame)],
            "tumor_cells": [int(labels.sum())],
            "normal_cells": [int((~labels).sum())],
            "candidate_partners": [len(partner_table)],
            "eligible_partners": [len(eligible_genes)],
            "train_patients": [len(np.unique(patients[train]))],
            "test_patients": [len(np.unique(patients[test]))],
        }
    ).to_csv(args.output_dir / "run_summary.csv", index=False)
    print(external.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
