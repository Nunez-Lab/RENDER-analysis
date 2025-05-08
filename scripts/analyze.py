################################################################################
# %% Imports

import matplotlib.pyplot as plt
import polars as pl

import glob
import os
import sys

################################################################################
# %% Command-line arguments

METADATA_DIR = sys.argv[1]
RNA_COUNTS_PATH = sys.argv[2]
RNA_DGE_DIR = sys.argv[3]
OUTPUT_DIR = sys.argv[4]

################################################################################
# %% RNA-seq count plots

# %% Load data

METADATA_DIR = "metadata"
RNA_COUNTS_PATH = "output/RNAseq/aggregated-reads/counts.csv"

cd55_prefix = "ENSG00000196352"
counts = pl.read_csv(RNA_COUNTS_PATH)

# %% Create RNA-seq count plots

rmeta = pl.read_csv(os.path.join(METADATA_DIR, "rna.csv"))

for (cell_line, condition), g in rmeta.group_by("cell_line", "condition"):
    assert len(g) == 2, "expected 2 replicates"

    rep1 = g.row(by_predicate=pl.col("replicate") == 1, named=True)
    rep2 = g.row(by_predicate=pl.col("replicate") == 2, named=True)

    print(rep1)

