################################################################################
# %% Imports

import importlib
import os
import sys

import matplotlib.ticker as mtick
import numpy as np
import polars as pl

try:
    import lib
except ModuleNotFoundError:
    import scripts.lib as lib

importlib.reload(lib)

################################################################################
# %% Command-line arguments

METADATA_DIR = sys.argv[1]
GENE_METADATA_DIR = sys.argv[2]

RNA_ABUNDANCE_PATH = sys.argv[3]
RNA_TEST_DIR = sys.argv[4]

EM_AVG_PATH = sys.argv[5]
EM_TEST_DIR = sys.argv[6]
EM_AGG_DIR = sys.argv[7]

OUTPUT_DIR = sys.argv[8]

################################################################################
# %% RNA-seq count plots

METADATA_DIR = "metadata"
GENE_METADATA_DIR = "output/gene-metadata"

RNA_ABUNDANCE_PATH = "output/RNAseq/aggregated-reads/abundance.csv"
RNA_TEST_DIR = "output/RNAseq/deseq2"

EM_AVG_PATH = "output/EMseq/average-methylation-info/info.tsv"
EM_TEST_DIR = "output/EMseq/dss"
EM_AGG_DIR = "output/EMseq/dss-aggregated"

OUTPUT_DIR = "output/analysis"

CD55_ENSG = "ENSG00000196352"
CLTA_ENSG = "ENSG00000122705"

# %% Load data

gene_metadata = pl.read_csv(
    os.path.join(GENE_METADATA_DIR, "gene-metadata.tsv"),
    separator="\t",
)

comparisons = pl.read_csv(os.path.join(METADATA_DIR, "comparisons.csv"))
rmeta = pl.read_csv(os.path.join(METADATA_DIR, "rna.csv"))
rabundance = pl.read_csv(RNA_ABUNDANCE_PATH).rename({"": "ensembl_gene_id_version"})

highlights = {
    "jurkat": pl.col("ensembl_gene_id_version").str.starts_with(CD55_ENSG),
    "hek293t": pl.col("ensembl_gene_id_version").str.starts_with(CLTA_ENSG),
}

# %% Create RNA-seq count plots for replicates

importlib.reload(lib)

for (cell_line, condition), g in rmeta.group_by("cell_line", "condition"):
    assert len(g) == 2, "expected 2 replicates"

    rep1 = g.row(by_predicate=pl.col("replicate") == 1, named=True)["sample"]
    rep2 = g.row(by_predicate=pl.col("replicate") == 2, named=True)["sample"]

    lib.rna_count_plot(
        rabundance,
        rep1,
        rep2,
        title=f"Replicates for {cell_line} {condition}",
    )[0].save(
        os.path.join(
            OUTPUT_DIR,
            f"rna-replicates-{cell_line}-{condition}.pdf",
        )
    )

# %% Create RNA-seq count plots for comparisons

importlib.reload(lib)

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]

    control_samples = rmeta.filter(
        pl.col("cell_line") == cell_line,
        pl.col("condition") == control,
    )["sample"]

    treatment_samples = rmeta.filter(
        pl.col("cell_line") == cell_line,
        pl.col("condition") == treatment,
    )["sample"]

    df = rabundance.select(
        pl.col("ensembl_gene_id_version"),
        pl.mean_horizontal(pl.col(control_samples)).alias("control"),
        pl.mean_horizontal(pl.col(treatment_samples)).alias("treatment"),
    )

    lib.rna_count_plot(
        df,
        "control",
        "treatment",
        title=f"RENDER {treatment} vs {control} for {cell_line}",
        highlight=highlights[cell_line],
    )[0].save(
        os.path.join(
            OUTPUT_DIR,
            f"rna-comparison-counts-{cell_line}-{control}-{treatment}.pdf",
        )
    )


# %% Create RNA-seq volcano plots for comparisons

importlib.reload(lib)

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]

    rtest = (
        pl.read_csv(
            os.path.join(
                RNA_TEST_DIR,
                f"{cell_line}-{control}-{treatment}.csv",
            ),
            null_values="NA",
        )
        .rename({"": "ensembl_gene_id_version"})
        .join(
            gene_metadata[["ensembl_gene_id_version", "external_gene_name"]],
            on="ensembl_gene_id_version",
            how="left",
            validate="1:1",
        )
    )

    lib.volcano_plot(
        rtest,
        title=f"RENDER {treatment} vs {control} for {cell_line}",
        highlight=highlights[cell_line],
        threshold=15,
        gene_name_feature="external_gene_name",
    )[0].save(
        os.path.join(
            OUTPUT_DIR,
            f"rna-volcano-{cell_line}-{control}-{treatment}.pdf",
        )
    )

# %% Load methylation data

importlib.reload(lib)

etest = {}

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]
    base = f"{cell_line}-{control}-{treatment}"

    etest[base] = lib.load_dss_results(
        os.path.join(EM_TEST_DIR, f"{base}.csv"),
    )

# % % Create EM-seq control plots

importlib.reload(lib)

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]
    base = f"{cell_line}-{control}-{treatment}"

    for condition in ["control", "treatment"]:
        feature = f"mu_{condition}"
        lib.manhattan(
            etest[base].filter(
                pl.col("chr_order") > 26,
            ),
            by="nice_chr",
            feature=feature,
            two_sided=None,
            highlight=pl.col(feature) > 0.5,
            highlight_color=lib.PURPLE,
            yticks=np.arange(0, 1.1, 0.2),
            ylabel=r"$\bf{\%}$ $\bf{reads}$ $\bf{methylated}$",
            yaxis_formatter=mtick.PercentFormatter(1.0),
            xtick_rotation=0,
        )[0].save(
            os.path.join(
                OUTPUT_DIR,
                f"em-control-{base}-{condition}.png",
            ),
            dpi=300,
        )

# % % Create main Manhattan plots

importlib.reload(lib)

score_yticks = np.arange(-250, 251, 50)
score_yticklabels = abs(score_yticks)

targeted_genes = {
    "jurkat": "CD55",
    "hek293t": "CLTA",
}

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]
    base = f"{cell_line}-{control}-{treatment}"

    targeted_gene = targeted_genes[cell_line]
    tss = lib.gene_info(gene_metadata, etest[base], targeted_gene)["tss_index"]

    print(f"Working on {base} Manhattan plot...")

    lib.manhattan(
        etest[base].filter(pl.col("chr_order") < 25),
        by="chr",
        feature="score",
        two_sided=(treatment, control),
        highlight=None,
        yticks=score_yticks,
        yticklabels=score_yticklabels,
        ylabel=r"$\bf{Significance}$ $\bf{score}$" + "\n" + r"$-\log_{10}($FDR$)$",
        use_xticks=True,
        arrow=(targeted_gene, tss),
    )[0].save(
        os.path.join(
            OUTPUT_DIR,
            f"em-manhattan-{base}.png",
        ),
        dpi=300,
    )

# %% Create correlation plots


importlib.reload(lib)

for row in comparisons.iter_rows(named=True):
    cell_line = row["cell_line"]
    control = row["control"]
    treatment = row["treatment"]
    base = f"{cell_line}-{control}-{treatment}"

    data = (
        pl.read_csv(
            os.path.join(
                RNA_TEST_DIR,
                f"{cell_line}-{control}-{treatment}.csv",
            ),
            null_values="NA",
        )
        .select(
            pl.col("").alias("ensembl_gene_id_version"),
            pl.col("padj").log(base=10).neg().alias("rscore"),
        )
        .join(
            pl.read_csv(
                os.path.join(
                    EM_AGG_DIR,
                    f"{cell_line}-{control}-{treatment}.csv",
                ),
            ).select(
                pl.col("ensembl_gene_id_version"),
                pl.col("min_fdr").log(base=10).neg().alias("escore"),
            ),
            on="ensembl_gene_id_version",
            how="inner",
            validate="1:1",
        )
        .join(
            gene_metadata[["ensembl_gene_id_version", "external_gene_name"]],
            on="ensembl_gene_id_version",
            how="left",
            validate="1:1",
        )
    )

    lib.correlation_plot(
        data,
        "rscore",
        "escore",
        xlabel=r"$\bf{RNAseq}$ $\bf{significance}$ $\bf{score}$"
        + "\n"
        + r"$-\log_{10}($adjusted $p$-value$)$",
        ylabel=r"$\bf{EMseq}$ $\bf{significance}$ $\bf{score}$"
        + "\n"
        + r"$-\log_{10}($FDR$)$",
        highlight=highlights[cell_line],
    )[0].save(
        os.path.join(
            OUTPUT_DIR,
            f"correlation-{base}.png",
        ),
        dpi=300,
    )
