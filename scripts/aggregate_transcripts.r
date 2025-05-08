################################################################################
# %% Imports

library(org.Hs.eg.db)
library(biomaRt)

################################################################################
# %% Command-line arguments

args = commandArgs(trailingOnly = TRUE)
METADATA_FILE = args[1]
KALLISTO_QUANT_DIR = args[2]
OUTPUT_DIR = args[3]

################################################################################
# %% Main script

# Set up file paths

metadata = read.csv(METADATA_FILE, header = TRUE)

files = file.path(
    KALLISTO_QUANT_DIR,
    metadata$sample,
    "abundance.h5"
)

names(files) = metadata$sample

# %% Create tx2gene using biomaRt

mart = useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 114,
)

tx2gene = getBM(
    attributes =
        c("ensembl_transcript_id_version",
          "ensembl_gene_id_version"),
    mart = mart,
)

# %% Aggregate abundances and save count file

txi = tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene
)

write.csv(txi$counts, OUTPUT_PATH)
