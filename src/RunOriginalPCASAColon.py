"""Run the original PCASA gate evaluator on downloaded GSE178341 CRC cells.

The upstream R random-forest stage is replaced with sklearn only to produce
the ranked 100-gene file consumed by the unchanged PCASA step-2 script.
Both the legacy outputs and corrected coverages are retained for comparison.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

from ModernPCASA import candidate_genes
from ValidateColonSingleCell import read_hacohen_target_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pcasa-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=Path("data/single_cell"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/pcasa_original_colon")
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--rf-cells", type=int, default=50_000)
    return parser.parse_args()


def corrected_gates(
    detected: np.ndarray, labels: np.ndarray, genes: list[str]
) -> pd.DataFrame:
    tumor, normal = detected[labels], detected[~labels]
    rows: list[dict[str, object]] = []
    for logic, gates in (
        ("AND", combinations(range(len(genes)), 2)),
        ("OR", combinations(range(len(genes)), 2)),
        ("NOT", permutations(range(len(genes)), 2)),
    ):
        for left, right in gates:
            if logic == "AND":
                pt, pn = tumor[:, left] & tumor[:, right], normal[:, left] & normal[:, right]
            elif logic == "OR":
                pt, pn = tumor[:, left] | tumor[:, right], normal[:, left] | normal[:, right]
            else:
                pt, pn = tumor[:, left] & ~tumor[:, right], normal[:, left] & ~normal[:, right]
            rows.append(
                {
                    "logic": logic,
                    "gate": f"{genes[left]}_{genes[right]}",
                    "tumor_coverage": pt.mean(),
                    "normal_coverage": pn.mean(),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_code = args.output_dir / "code"
    run_data = args.output_dir / "data"
    run_results = args.output_dir / "results"
    for directory in (run_code, run_data, run_results):
        directory.mkdir(exist_ok=True)

    genes = candidate_genes(Path("data/test-normalized-matrix.txt"))
    stem = "GSE178341_crc10x_full_c295v4_submit"
    counts, barcodes = read_hacohen_target_counts(
        args.data_dir / f"{stem}.h5", genes, 10_000
    )
    clusters = pd.read_csv(args.data_dir / f"{stem}_cluster.csv.gz")
    clinical = pd.read_csv(args.data_dir / f"{stem}_metatables.csv.gz")
    if not np.array_equal(barcodes, clusters["sampleID"].to_numpy()):
        raise ValueError("Expression and annotation barcodes are not aligned")
    keep = clusters["clTopLevel"].eq("Epi").to_numpy() & clinical[
        "SPECIMEN_TYPE"
    ].isin(["T", "N"]).to_numpy()
    expression = counts.to_numpy()[keep]
    detected = expression > 0
    labels = clinical.loc[keep, "SPECIMEN_TYPE"].eq("T").to_numpy()
    kept_barcodes = barcodes[keep].astype(str)

    rng = np.random.default_rng(args.seed)
    sampled = np.concatenate(
        [
            rng.choice(np.flatnonzero(labels), min(args.rf_cells // 2, labels.sum()), False),
            rng.choice(
                np.flatnonzero(~labels),
                min(args.rf_cells // 2, (~labels).sum()),
                False,
            ),
        ]
    )
    forest = RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced_subsample",
        n_jobs=-1,
        random_state=args.seed,
    )
    forest.fit(detected[sampled], labels[sampled])
    order = np.argsort(forest.feature_importances_)[::-1][:100]
    ranked_genes = [genes[index] for index in order]

    annotation = pd.DataFrame(
        {
            "CellCode": kept_barcodes,
            "CellType": np.where(labels, "Tumor", "Normal"),
        }
    )
    annotation.to_csv(run_data / "input-2a__scrna_annotation.txt", sep="\t", index=False)
    matrix = pd.DataFrame(
        detected[:, order].T.astype(np.uint8),
        index=ranked_genes,
        columns=kept_barcodes,
    )
    matrix.index.name = "Gene_symbol"
    matrix.to_csv(run_data / "input-2b__scrna_gc-matrix.txt", sep="\t")
    with (run_results / "result-1__rf_var-imp.txt").open("w", newline="\n") as handle:
        handle.write('"FALSE" "TRUE" "MeanDecreaseAccuracy" "MeanDecreaseGini"\n')
        for rank, index in enumerate(order, 1):
            importance = forest.feature_importances_[index]
            handle.write(
                f'"{genes[index]}" 0 0 {importance:.12g} {100-rank:.1f}\n'
            )

    original = args.pcasa_dir / "code" / "step-2__gate_coverage_calc.py"
    copied = run_code / original.name
    shutil.copy2(original, copied)
    completed = subprocess.run(
        [
            sys.executable,
            copied.name,
            str((run_data / "input-2a__scrna_annotation.txt").resolve()),
            str((run_data / "input-2b__scrna_gc-matrix.txt").resolve()),
        ],
        cwd=run_code,
        check=True,
        capture_output=True,
        text=True,
    )
    (args.output_dir / "legacy_step2.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )

    corrected = corrected_gates(detected[:, order], labels, ranked_genes)
    corrected.to_csv(run_results / "corrected_gate_coverage.csv", index=False)
    pd.DataFrame(
        {
            "gene": ranked_genes,
            "rf_importance": forest.feature_importances_[order],
        }
    ).to_csv(run_results / "python_rf_ranking.csv", index=False)
    print(
        f"Tested {len(labels):,} cells: {labels.sum():,} tumor and "
        f"{(~labels).sum():,} normal; outputs: {run_results}"
    )


if __name__ == "__main__":
    main()
