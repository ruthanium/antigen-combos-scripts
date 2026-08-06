"""Python port of CalculateClusteringScores.R."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pipeline_core import (
    add_antigen_labels,
    all_combinations,
    cancer_normal_subset,
    cancer_types,
    clustering_scores,
    gene_columns,
    read_expression,
    read_sketch,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("data/test-normalized-matrix.txt"))
    parser.add_argument("--clinical", type=Path, default=Path("data/clin-antigens.txt"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--test-run", action="store_true", help="Limit doubles/triples to 50 per cancer")
    parser.add_argument("--cancer", action="append", help="Only process this cancer (repeatable)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    matrix = read_expression(args.matrix)
    clinical = set(pd.read_csv(args.clinical, header=None)[0].astype(str))
    genes = gene_columns(matrix)
    cancers = args.cancer or cancer_types(matrix)
    limit = 50 if args.test_run else None

    single_frames = []
    double_frames = []
    for cancer in cancers:
        print(f"{cancer}: singles/doubles", flush=True)
        sketch = read_sketch(args.results_dir, cancer)
        subset = cancer_normal_subset(matrix, cancer, sketch)
        single_frames.append(clustering_scores(subset, cancer, ((g,) for g in genes)))
        double_frames.append(
            clustering_scores(subset, cancer, all_combinations(genes, 2, limit))
        )

    singles = add_antigen_labels(pd.concat(single_frames, ignore_index=True), clinical, 1)
    doubles = add_antigen_labels(pd.concat(double_frames, ignore_index=True), clinical, 2)
    singles.to_csv(args.results_dir / "singles.csv", index=False)
    doubles.to_csv(args.results_dir / "doubles.csv", index=False)

    filtered = singles.loc[(singles["db"] <= 5) & (singles["dist.man"] > 2)]
    triple_frames = []
    for cancer in cancers:
        candidate_genes = filtered.loc[filtered["can"].eq(cancer), "combo"].tolist()
        print(f"{cancer}: triples from {len(candidate_genes)} genes", flush=True)
        sketch = read_sketch(args.results_dir, cancer)
        subset = cancer_normal_subset(matrix, cancer, sketch)
        triple_frames.append(
            clustering_scores(subset, cancer, all_combinations(candidate_genes, 3, limit))
        )

    nonempty = [frame for frame in triple_frames if not frame.empty]
    triples = pd.concat(nonempty, ignore_index=True) if nonempty else pd.DataFrame(
        columns=["can", "combo", "db", "dist.man", "ttype"]
    )
    triples = add_antigen_labels(triples, clinical, 3)
    triples.to_csv(args.results_dir / "triples.csv", index=False)
    print(
        f"Wrote {len(singles)} singles, {len(doubles)} doubles, and {len(triples)} triples",
        flush=True,
    )


if __name__ == "__main__":
    main()
