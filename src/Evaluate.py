"""Python port of Evaluate.R."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pipeline_core import evaluate_combo, read_expression, read_sketch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("data/test-normalized-matrix.txt"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    return parser.parse_args()


def evaluate_table(
    matrix: pd.DataFrame, scores: pd.DataFrame, results_dir: Path, size: int
) -> pd.DataFrame:
    if size == 1:
        scores = scores.loc[scores["ttype"].eq("H")].copy()
    else:
        scores = scores.loc[~scores["ttype"].eq("L" * size)].copy()

    sketches: dict[str, set[str]] = {}
    metrics = []
    for index, row in scores.iterrows():
        cancer = row["can"]
        if cancer not in sketches:
            sketches[cancer] = read_sketch(results_dir, cancer)
        result = evaluate_combo(matrix, cancer, row["combo"].split(":"), sketches[cancer])
        result["_index"] = index
        metrics.append(result)
    metric_frame = pd.DataFrame(metrics).set_index("_index") if metrics else pd.DataFrame()
    return scores.join(metric_frame)


def main() -> None:
    args = parse_args()
    matrix = read_expression(args.matrix)
    for name, size in (("singles", 1), ("doubles", 2), ("triples", 3)):
        source = args.results_dir / f"{name}.csv"
        scores = pd.read_csv(source)
        print(f"Evaluating {len(scores)} {name}", flush=True)
        total = evaluate_table(matrix, scores, args.results_dir, size)
        total.to_csv(args.results_dir / f"{name}.total.csv", index=False)
        print(f"Wrote {len(total)} rows to {name}.total.csv", flush=True)


if __name__ == "__main__":
    main()
