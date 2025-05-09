################################################################################
# %% Imports

library(DESeq2)

################################################################################
# %% Command-line arguments

args = commandArgs(trailingOnly = TRUE)
GENE_METADATA_FILE = args[1]
METADATA_FILE = args[2]
COMPARISONS_FILE = args[3]
COUNTS_PATH = args[4]
OUTPUT_PATH = args[5]

################################################################################
# %% Main script

gene_metadata = read.table(
    GENE_METADATA_FILE,
    header=TRUE,
    sep="\t",
)

protein_coding = gene_metadata[
    gene_metadata$gene_biotype == "protein_coding",
    "ensembl_gene_id_version"
]

metadata = read.csv(
    METADATA_FILE,
    header=TRUE,
)

comparisons = read.csv(
    COMPARISONS_FILE,
    header=TRUE,
)

counts = round(
    read.csv(
        COUNTS_PATH,
        row.names=1
    )
)

counts = counts[
    rownames(counts) %in% protein_coding,
]

# %%

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
        condition = as.factor(conditions)
    )

    # Run DESeq2

    dds = DESeqDataSetFromMatrix(
        countData = counts[, cols],
        colData = colData,
        design = ~ condition
    )

    dds = DESeq(dds)

    res = results(
        dds,
        contrast=c("condition", row$treatment, row$control)
    )

    # Write results to CSV

    filename = paste0(
        paste(row$cell_line, row$control, row$treatment, sep="-"),
        ".csv"
    )

    write.csv(
        res,
        file.path(OUTPUT_PATH, filename),
    )
}
