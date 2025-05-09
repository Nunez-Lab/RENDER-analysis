################################################################################
# %% Imports

import glob
import importlib
import os
import sys

import matplotlib.ticker as mtick
import numpy as np
import polars as pl

################################################################################
# %% Command-line arguments

GENE_METADATA_PATH = sys.argv[1]
DSS_DIR = sys.argv[2]
OUTPUT_DIR = sys.argv[3]

################################################################################
# %% Main script

GENE_METADATA_PATH = "output/gene-metadata/gene-metadata.tsv"
DSS_DIR = "output/EMseq/dss"
OUTPUT_DIR = "output/EMseq/dss-aggregated"

gene_metadata = (
    pl.scan_csv(GENE_METADATA_PATH, separator="\t")
    .with_columns(
        promoter_start=pl.col("start_position") - 500,
        promoter_end=pl.col("end_position") + 500,
    )
    .select(
        "ensembl_gene_id_version",
        "chromosome_name",
        "promoter_start",
        "promoter_end",
    )
)

for methyl_path in sorted(glob.glob(f"{DSS_DIR}/*.csv")):
    (
        pl.scan_csv(methyl_path)
        .filter(~pl.col("chr").str.contains("_"))
        .with_columns(
            chromosome_name=pl.when(pl.col("chr") == "chrX")
            .then(pl.lit("X"))
            .when(pl.col("chr") == "chrY")
            .then(pl.lit("Y"))
            .when(pl.col("chr") == "chrM")
            .then(pl.lit("MT"))
            .otherwise(pl.col("chr").str.slice(3)),
        )
        .select(
            "chromosome_name",
            "pos",
            "pval",
            "fdr",
            "diff",
        )
        .sort(by="pos")
        .join_asof(
            gene_metadata.sort(by="promoter_start"),
            by="chromosome_name",
            left_on="pos",
            right_on="promoter_start",
            strategy="backward",
            coalesce=False,
        )
        .filter(~pl.col("promoter_end").is_null())
        .filter(pl.col("pos") <= pl.col("promoter_end"))
        .group_by("ensembl_gene_id_version")
        .agg(
            min_pval=pl.col("pval").min(),
            min_fdr=pl.col("fdr").min(),
            mean_diff=pl.col("diff").mean(),
        )
        .collect()
        .write_csv(
            os.path.join(
                OUTPUT_DIR,
                os.path.basename(methyl_path),
            ),
        )
    )
