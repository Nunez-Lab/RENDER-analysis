import polars as pl
import matplotlib.pyplot as plt


import os
import sys
import glob

RESULTS_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]


CD55_IDS = [
    "NM_000574.5",
    "NM_001114752.3",
    "NM_001300902.2",
    "NM_001300903.2",
    "NM_001300904.2",
]

for results_path in sorted(glob.glob(f"{RESULTS_DIR}/*.csv")):
    plot_name = os.path.basename(results_path)

    # Load in the abundance data table (separated by tabs, not commas)
    results = pl.read_csv(
        results_path,
        null_values=["NA"],
    ).rename(
        {"": "target_id"}
    ).with_columns(
        score=-pl.col("padj").log(base=10),
        color=pl.when(pl.col("target_id").is_in(CD55_IDS)).then(pl.lit("red")).otherwise(pl.lit("gray"))
    ).filter(pl.col("target_id").str.starts_with("NM_"))

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    for (c,), g in results.group_by("color"):
        if c == "red":
            alpha = 1
            zorder = 10
        else:
            alpha = 0.1
            zorder = 5
        ax.scatter(g["log2FoldChange"], g["score"], color=c, alpha=alpha, zorder=zorder)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR + "/" + plot_name + ".png")

    # print(results.top_k(10, by=pl.col("score")))
