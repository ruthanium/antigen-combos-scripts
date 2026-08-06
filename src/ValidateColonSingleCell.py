"""Validate colon-cancer antigen candidates in two public single-cell atlases.

Sources
-------
1. Pelka/Hacohen CRC atlas, GEO GSE178341.
2. ICBI integrated CRC atlas, https://crc.icbi.at/.

The script reports detection fractions and mean expression by cell compartment.
It deliberately distinguishes "tumor epithelial" (a specimen/lineage label in
GSE178341) from the explicit "Cancer cell" annotation in the ICBI atlas.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy.sparse import csc_matrix
from server.common.fbs.matrix import decode_matrix_fbs


DEFAULT_GENES = ("GRIN2D", "RNF43", "NOX1", "ABCA12", "AGTR1")
DEFAULT_COMBOS = (("ABCA12", "NOX1"), ("AGTR1", "GRIN2D", "NOX1"))
ICBI_API = "https://crc.icbi.at/atlas/api/v0.2/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data/single_cell"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/single_cell_validation")
    )
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--h5-chunk-cells", type=int, default=10_000)
    return parser.parse_args()


def summarize(
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    group_columns: list[str],
    patient_column: str,
    source: str,
) -> pd.DataFrame:
    rows = []
    joined = pd.concat([metadata.reset_index(drop=True), expression], axis=1)
    for group_values, group in joined.groupby(group_columns, observed=True, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        group_info = dict(zip(group_columns, group_values))
        for gene in expression.columns:
            values = group[gene].to_numpy(dtype=float)
            positive = values > 0
            patient_rates = (
                group.assign(_positive=positive)
                .groupby(patient_column, observed=True)["_positive"]
                .mean()
            )
            rows.append(
                {
                    "source": source,
                    **group_info,
                    "gene": gene,
                    "n_cells": len(values),
                    "pct_expressed": 100 * float(positive.mean()),
                    "mean_expression": float(values.mean()),
                    "mean_when_expressed": (
                        float(values[positive].mean()) if positive.any() else 0.0
                    ),
                    "n_patients": int(patient_rates.size),
                    "median_patient_pct_expressed": (
                        100 * float(patient_rates.median()) if patient_rates.size else np.nan
                    ),
                    "pct_patients_with_expression": (
                        100 * float((patient_rates > 0).mean())
                        if patient_rates.size
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def summarize_combinations(
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    group_columns: list[str],
    patient_column: str,
    source: str,
) -> pd.DataFrame:
    rows = []
    joined = pd.concat([metadata.reset_index(drop=True), expression], axis=1)
    for group_values, group in joined.groupby(group_columns, observed=True, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        group_info = dict(zip(group_columns, group_values))
        for genes in DEFAULT_COMBOS:
            if not set(genes).issubset(expression.columns):
                continue
            positive = group.loc[:, list(genes)].gt(0).all(axis=1)
            patient_rates = (
                group.assign(_positive=positive)
                .groupby(patient_column, observed=True)["_positive"]
                .mean()
            )
            rows.append(
                {
                    "source": source,
                    **group_info,
                    "combo": "+".join(genes),
                    "n_cells": len(group),
                    "n_combo_positive": int(positive.sum()),
                    "pct_combo_positive": 100 * float(positive.mean()),
                    "n_patients": int(patient_rates.size),
                    "median_patient_pct_combo_positive": (
                        100 * float(patient_rates.median()) if patient_rates.size else np.nan
                    ),
                    "pct_patients_with_combo": (
                        100 * float((patient_rates > 0).mean())
                        if patient_rates.size
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def read_hacohen_target_counts(
    h5_path: Path, genes: list[str], chunk_cells: int
) -> tuple[pd.DataFrame, np.ndarray]:
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["matrix"]
        feature_names = matrix["features/name"][:].astype(str)
        feature_lookup = {name: index for index, name in enumerate(feature_names)}
        missing = [gene for gene in genes if gene not in feature_lookup]
        if missing:
            raise ValueError(f"Genes absent from GSE178341: {missing}")
        target_rows = np.array([feature_lookup[gene] for gene in genes])
        shape = tuple(int(value) for value in matrix["shape"][:])
        n_features, n_cells = shape
        indptr = matrix["indptr"][:]
        # GEO stored each sparse array as one very large gzip-compressed HDF5
        # chunk. Reading it in cell slices would decompress the same 9 GB of
        # arrays repeatedly, so load each array exactly once and immediately
        # subset the resulting CSC matrix to the requested genes.
        print("GSE178341 counts: decompressing sparse arrays once", flush=True)
        sparse = csc_matrix(
            (matrix["data"][:], matrix["indices"][:], indptr),
            shape=(n_features, n_cells),
        )
        output = sparse[target_rows, :].T.toarray().astype(np.float32, copy=False)
        del sparse
        print(f"GSE178341 counts: extracted {n_cells:,} cells", flush=True)

        barcodes = matrix["barcodes"][:].astype(str)
    return pd.DataFrame(output, columns=genes), barcodes


def analyze_hacohen(data_dir: Path, output_dir: Path, genes: list[str], chunk: int):
    stem = "GSE178341_crc10x_full_c295v4_submit"
    counts, barcodes = read_hacohen_target_counts(
        data_dir / f"{stem}.h5", genes, chunk
    )
    clusters = pd.read_csv(data_dir / f"{stem}_cluster.csv.gz")
    clinical = pd.read_csv(data_dir / f"{stem}_metatables.csv.gz")
    if not np.array_equal(barcodes, clusters["sampleID"].to_numpy()):
        raise ValueError("GSE178341 H5 barcodes do not align with cluster annotations")
    if not np.array_equal(barcodes, clinical["cellID"].to_numpy()):
        raise ValueError("GSE178341 H5 barcodes do not align with clinical annotations")

    metadata = pd.DataFrame(
        {
            "specimen": clinical["SPECIMEN_TYPE"].map({"T": "tumor", "N": "normal"}),
            "cell_type": clusters["clTopLevel"],
            "cell_subtype": clusters["cl295v11SubFull"],
            "patient_id": clinical["PID"],
        }
    )
    summary = summarize(
        metadata,
        counts,
        ["specimen", "cell_type"],
        "patient_id",
        "GSE178341 Pelka/Hacohen",
    )
    summary.to_csv(output_dir / "hacohen_expression_by_compartment.csv", index=False)
    combo_summary = summarize_combinations(
        metadata,
        counts,
        ["specimen", "cell_type"],
        "patient_id",
        "GSE178341 Pelka/Hacohen",
    )
    combo_summary.to_csv(
        output_dir / "hacohen_combinations_by_compartment.csv", index=False
    )
    return summary


def get_icbi_metadata() -> pd.DataFrame:
    annotations = (
        "annotation-name=atlas_cell_type_coarse"
        "&annotation-name=sample_type"
        "&annotation-name=patient_id"
        "&annotation-name=dataset"
    )
    response = requests.get(
        ICBI_API + "annotations/obs?" + annotations,
        headers={"Accept": "application/octet-stream"},
        timeout=240,
    )
    response.raise_for_status()
    return decode_matrix_fbs(response.content).reset_index(drop=True)


def get_icbi_gene(gene: str) -> np.ndarray:
    response = requests.get(
        ICBI_API + f"data/var?var:var_names={gene}",
        headers={"Accept": "application/octet-stream"},
        timeout=240,
    )
    response.raise_for_status()
    frame = decode_matrix_fbs(response.content)
    if frame.shape[1] != 1:
        raise ValueError(f"Expected one ICBI expression column for {gene}")
    return frame.iloc[:, 0].to_numpy(dtype=np.float32)


def analyze_icbi(output_dir: Path, genes: list[str]):
    metadata = get_icbi_metadata()
    expression = pd.DataFrame(index=metadata.index)
    for gene in genes:
        print(f"ICBI atlas: fetching {gene}", flush=True)
        expression[gene] = get_icbi_gene(gene)
    if len(metadata) != len(expression):
        raise ValueError("ICBI metadata and expression row counts do not match")
    summary = summarize(
        metadata,
        expression,
        ["sample_type", "atlas_cell_type_coarse"],
        "patient_id",
        "ICBI integrated CRC atlas",
    )
    summary.to_csv(output_dir / "icbi_expression_by_compartment.csv", index=False)
    combo_summary = summarize_combinations(
        metadata,
        expression,
        ["sample_type", "atlas_cell_type_coarse"],
        "patient_id",
        "ICBI integrated CRC atlas",
    )
    combo_summary.to_csv(
        output_dir / "icbi_combinations_by_compartment.csv", index=False
    )
    return summary


def make_plot(hacohen: pd.DataFrame, icbi: pd.DataFrame, output_dir: Path) -> None:
    hacohen_focus = hacohen.loc[
        ((hacohen["specimen"] == "tumor") & hacohen["cell_type"].isin(["Epi", "Strom", "Myeloid", "TNKILC"]))
        | ((hacohen["specimen"] == "normal") & hacohen["cell_type"].eq("Epi"))
    ].copy()
    hacohen_focus["compartment"] = (
        hacohen_focus["specimen"] + " " + hacohen_focus["cell_type"]
    )

    icbi_focus = icbi.loc[
        (
            icbi["sample_type"].eq("primary tumor")
            & icbi["atlas_cell_type_coarse"].isin(
                ["Cancer cell", "Epithelial cell", "Stromal cell", "Myeloid cell", "T cell"]
            )
        )
        | (
            icbi["sample_type"].isin(["adjacent normal", "healthy normal"])
            & icbi["atlas_cell_type_coarse"].eq("Epithelial cell")
        )
    ].copy()
    icbi_focus["compartment"] = (
        icbi_focus["sample_type"] + " " + icbi_focus["atlas_cell_type_coarse"]
    )

    plot_data = pd.concat(
        [
            hacohen_focus[["source", "compartment", "gene", "pct_expressed", "mean_expression"]],
            icbi_focus[["source", "compartment", "gene", "pct_expressed", "mean_expression"]],
        ],
        ignore_index=True,
    )
    plot_data.to_csv(output_dir / "focused_compartment_summary.csv", index=False)

    sources = plot_data["source"].drop_duplicates().tolist()
    fig, axes = plt.subplots(1, len(sources), figsize=(9 * len(sources), 5), squeeze=False)
    for axis, source in zip(axes[0], sources):
        subset = plot_data.loc[plot_data["source"].eq(source)]
        order = subset["compartment"].drop_duplicates().tolist()
        sns.scatterplot(
            data=subset,
            x="gene",
            y="compartment",
            size="pct_expressed",
            hue="mean_expression",
            sizes=(15, 350),
            palette="viridis",
            edgecolor="black",
            linewidth=0.25,
            ax=axis,
        )
        axis.set_title(source)
        axis.set(xlabel="", ylabel="")
        axis.set_yticks(range(len(order)), order)
        axis.tick_params(axis="x", rotation=45)
        axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "colon_candidate_single_cell_dotplot.png", dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    hacohen = analyze_hacohen(
        args.data_dir, args.output_dir, args.genes, args.h5_chunk_cells
    )
    icbi = analyze_icbi(args.output_dir, args.genes)
    make_plot(hacohen, icbi, args.output_dir)
    print(f"Results written to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
