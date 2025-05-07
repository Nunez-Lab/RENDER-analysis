################################################################################
# Imports

library(DESeq2)

################################################################################
# Command-line arguments

args = commandArgs(trailingOnly = TRUE)
CTS_PATH = args[1]
OUTPUT_PATH = args[2]

################################################################################
# Main script

# TODO: Hard-coded samples of interest
cols = c("DX5N1_S8", "DX5N2_S21", "DX5S1_S9", "DX5S2_S22")
colData = data.frame(id = cols, condition = c("non_targeting", "non_targeting", "single", "single"))

# Load read count matrix
cts = as.matrix(read.csv(CTS_PATH, sep="\t", row.names="target_id"))[, cols]

# Run DESeq2
dds = DESeqDataSetFromMatrix(
    countData = cts,
    colData = colData,
    design = ~ condition
)
dds = DESeq(dds)
res = results(dds)

# Save results to a TSV
write.csv(res, OUTPUT_PATH)
