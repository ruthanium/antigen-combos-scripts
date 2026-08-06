"""Non-interactive Python port of the ranking portion of Analyze.R."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pipeline_core import clean_and_scale, diverse_top_k


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def normalized_label(label: str, size: int) -> str:
    parts = label.split(":")
    if size > 1:
        parts.sort()
    return ":".join(parts)


def top_table(frame: pd.DataFrame, size: int, k: int) -> pd.DataFrame:
    frame = clean_and_scale(frame)
    frame["lab"] = frame["lab"].map(lambda value: normalized_label(value, size))
    groups = []
    for (_, _), group in frame.groupby(["can", "lab"], sort=False):
        if size == 1:
            groups.append(group.head(k))
        else:
            groups.append(diverse_top_k(group, size=size, k=k))
    return pd.concat(groups, ignore_index=True) if groups else frame.iloc[0:0]


def main() -> None:
    args = parse_args()
    categories = []
    top_frames = []
    for name, size, category in (
        ("singles", 1, "single"),
        ("doubles", 2, "double"),
        ("triples", 3, "triple"),
    ):
        total = pd.read_csv(args.results_dir / f"{name}.total.csv")
        top = top_table(total, size, args.top_k)
        top["category"] = category
        top.to_csv(args.results_dir / f"top_{name}.csv", index=False)
        top_frames.append(top)
        categories.append(category)

    combined = pd.concat(top_frames, ignore_index=True)
    combined.to_csv(args.results_dir / "top_combinations.csv", index=False)

    if not combined.empty:
        sns.set_theme(style="whitegrid")
        fig, axis = plt.subplots(figsize=(7, 4))
        sns.boxplot(data=combined, x="category", y="f1", order=categories, ax=axis)
        axis.set(xlabel="", ylabel="F1 score", ylim=(0, 1))
        fig.tight_layout()
        fig.savefig(args.results_dir / "f1_by_gate_size.png", dpi=180)
        plt.close(fig)

        summary = (
            combined.groupby("category")[["f1", "prec", "rec"]]
            .agg(["count", "median", "mean"])
            .round(6)
        )
        summary.to_csv(args.results_dir / "performance_summary.csv")
    print(f"Wrote ranked tables for {len(combined)} top combinations", flush=True)


if __name__ == "__main__":
    main()
