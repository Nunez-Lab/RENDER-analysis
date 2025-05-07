import polars as pl

import glob
import os
import sys

QUANT_DIR = sys.argv[1]
RESULT_PATH = sys.argv[2]

dfs = []

for folder in sorted(glob.glob(f"{QUANT_DIR}/*")):
    sample_name = os.path.basename(folder)
    dfs.append(
        pl.read_csv(
            folder + "/abundance.tsv",
            separator="\t",
        ).select(
            pl.col("target_id"),
            pl.col("est_counts").round().cast(int).alias(sample_name)
        )
    )


df = dfs[0]

for df2 in dfs[1:]:
    df = df.join(df2, on="target_id", validate="1:1")

df.write_csv(RESULT_PATH, separator="\t")
