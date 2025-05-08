################################################################################
# %% Imports

library(DESeq2)

################################################################################
# %% Command-line arguments

args = commandArgs(trailingOnly = TRUE)
METADATA_FILE = args[1]
COMPARISONS_FILE = args[2]
COUNTS_PATH = args[3]
OUTPUT_PATH = args[4]

################################################################################
# %% Main script

metadata = read.csv(METADATA_FILE, header=TRUE)
comparisons = read.csv(COMPARISONS_FILE, header=TRUE)

for (i in 1:nrow(comparisons)) {
    # Set up DESeq2 input data

    row = comparisons[i,]
    control = metadata[
        (metadata$cell_line == row$cell_line) &
        (metadata$condition == row$control),
    ]
    treatment = metadata[
        (metadata$cell_line == row$cell_line) &
        (metadata$condition == row$treatment),
    ]

    cols = c(control$sample, treatment$sample)
    conditions = c(control$condition, treatment$condition)

    colData = data.frame(
        id = cols,
        condition = conditions
    )

    counts = as.matrix(
        read.csv(
            COUNTS_PATH,
            sep="\t",
            row.names="target_id"
        )
    )[, cols]

    # Run DESeq2

    dds = DESeqDataSetFromMatrix(
        countData = counts,
        colData = colData,
        design = ~ condition
    )
    dds = DESeq(dds)
    res = results(dds)

    # Write results to TSV

    filename = paste0(
        paste(row$cell_line, row$control, row$treatment, sep="-"),
        ".tsv"
    )

    write.table(
        res,
        file.path(OUTPUT_PATH, filename),
        quote=FALSE,
        sep="\t"
    )
}
