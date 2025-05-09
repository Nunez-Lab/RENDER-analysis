################################################################################
# %% Imports

library(biomaRt)
library(tidyverse)

################################################################################
# %% Command-line arguments

args = commandArgs(trailingOnly = TRUE)
OUTPUT_DIR = args[1]

################################################################################
# %% Main script

OUTPUT_DIR = "output/gene-metadata/"

mart = useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 114
)

tab = getBM(
    attributes = c(
          "ensembl_gene_id",
          "ensembl_gene_id_version",
          "chromosome_name",
          "start_position",
          "end_position",
          "external_gene_name"
    ),
    mart = mart
)

write_tsv(tab, file.path(OUTPUT_DIR, "gene-metadata.tsv"))
