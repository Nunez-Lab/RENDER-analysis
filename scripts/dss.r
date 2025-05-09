################################################################################
# %% Imports

library(tidyverse)
library(bsseq)
library(DSS)

################################################################################
# %% Command-line arguments

args = commandArgs(trailingOnly = TRUE)
METADATA_FILE = args[1]
COMPARISONS_FILE = args[2]
METHYLATION_DIR = args[3]
OUTPUT_PATH = args[4]

# %% Load metadata

metadata = read.csv(
    METADATA_FILE,
    header=TRUE,
)

comparisons = read.csv(
    COMPARISONS_FILE,
    header=TRUE,
)

# %% Load data

print("Loading data...")

data = list()

for (i in 1:nrow(metadata)) {
    sample = metadata[i, "sample"]
    data[[sample]] = as.data.frame(
        read_tsv(
            file.path(
                METHYLATION_DIR,
                paste0(sample, ".tsv")
            )
        )
    )
}

# %% Run differential methylation test

for (i in 1:nrow(comparisons)) {
    # Set up input data

    row = comparisons[i,]
    control = metadata[
        (metadata$cell_line == row$cell_line) &
        (metadata$condition == row$control),
    ]
    treatment = metadata[
        (metadata$cell_line == row$cell_line) &
        (metadata$condition == row$treatment),
    ]

    samples = c(control$sample, treatment$sample)
    conditions = c(control$condition, treatment$condition)

    dat = list()

    for (sample in samples) {
        dat[[sample]] = data[[sample]]
    }

    # Call DSS

    print("Making BSobj...")

    BSobj = makeBSseqData(
        dat,
        samples
    )

    print("Running differential methylation test...")

    dmlTest = DMLtest(
        BSobj,
        group1=control$sample,
        group2=treatment$sample,
        smoothing=TRUE,
        ncores=8
    )

    # Write results to CSV

    print("Writing CSV...")

    filename = paste0(
        paste(row$cell_line, row$control, row$treatment, sep="-"),
        ".csv"
    )

    write_csv(
        dmlTest,
        file.path(OUTPUT_PATH, filename),
    )
}
