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

# METADATA_DIR = "metadata"
# GENE_METADATA_DIR = "output/gene-metadata"
#
# RNA_ABUNDANCE_PATH = "output/RNAseq/aggregated-reads/abundance.csv"
# RNA_TEST_DIR = "output/RNAseq/deseq2"
#
# EM_AVG_PATH = "output/EMseq/average-methylation-info/info.tsv"
# EM_TEST_DIR = "output/EMseq/dss"
# EM_AGG_DIR = "output/EMseq/dss-aggregated"
#
# OUTPUT_DIR = "output/analysis"

# %% Load metadata

gene_metadata = pl.read_csv(
    os.path.join(GENE_METADATA_DIR, "gene-metadata.tsv"),
    separator="\t",
)

comparisons = pl.read_csv(os.path.join(METADATA_DIR, "comparisons.csv"))

targeted_genes = {
    "jurkat": "CD55",
    "hek293t": "CLTA",
}

nice_cell_lines = {
    "jurkat": "Jurkat",
    "hek293t": "HEK293T",
}


def comparisons_iter():
    for row in comparisons.iter_rows(named=True):
        cell_line = row["cell_line"]
        control = row["control"]
        treatment = row["treatment"]
        base = f"{cell_line}-{control}-{treatment}"
        yield cell_line, control, treatment, base


# %% Load RNA-seq data

rmeta = pl.read_csv(os.path.join(METADATA_DIR, "rna.csv"))

rabundance = (
    pl.read_csv(RNA_ABUNDANCE_PATH)
    .rename({"": "ensembl_gene_id_version"})
    .join(
        gene_metadata,
        on="ensembl_gene_id_version",
        how="left",
        validate="1:1",
    )
)

# %% Load methylation data

importlib.reload(lib)

etest = {}

for cell_line, control, treatment, base in comparisons_iter():
    etest[base] = lib.load_dss_results(
        os.path.join(EM_TEST_DIR, f"{base}.csv"),
    )

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
    )[0].save_organized(
        OUTPUT_DIR,
        "01-RNAseq-replicates",
        cell_line + "-" + condition,
    )

# %% Create RNA-seq count plots for comparisons

importlib.reload(lib)

for cell_line, control, treatment, base in comparisons_iter():
    control_samples = rmeta.filter(
        pl.col("cell_line") == cell_line,
        pl.col("condition") == control,
    )["sample"]

    treatment_samples = rmeta.filter(
        pl.col("cell_line") == cell_line,
        pl.col("condition") == treatment,
    )["sample"]

    df = rabundance.select(
        pl.col("external_gene_name"),
        pl.mean_horizontal(pl.col(control_samples)).alias("control"),
        pl.mean_horizontal(pl.col(treatment_samples)).alias("treatment"),
    )

    lib.rna_count_plot(
        df,
        "control",
        "treatment",
        xlabel=r"Non\!-\!targeting",
        ylabel=r"Targeting\ " + targeted_genes[cell_line],
        title="RENDER in " + nice_cell_lines[cell_line] + "s",
        highlight=pl.col("external_gene_name") == targeted_genes[cell_line],
    )[0].save_organized(
        OUTPUT_DIR,
        "02-RNAseq-comparison",
        base,
    )


# %% Create RNA-seq volcano plots for comparisons

importlib.reload(lib)

for cell_line, control, treatment, base in comparisons_iter():
    rtest = (
        pl.read_csv(
            os.path.join(RNA_TEST_DIR, f"{base}.csv"),
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
        treatment_name="targeting",
        control_name="non-targeting",
        highlight=pl.col("external_gene_name") == targeted_genes[cell_line],
        threshold=15,
        gene_name_feature="external_gene_name",
    )[0].save_organized(
        OUTPUT_DIR,
        "03-RNAseq-volcano",
        base,
    )

# %% Create EM-seq control plots

importlib.reload(lib)

for cell_line, control, treatment, base in comparisons_iter():
    for condition in ["control", "treatment"]:
        feature = f"mu_{condition}"
        lib.manhattan(
            etest[base].filter(pl.col("chr_order") > 26),
            by="nice_chr",
            feature=feature,
            two_sided=None,
            highlight=pl.col(feature) > 0.5,
            highlight_color=lib.PURPLE,
            yticks=np.arange(0, 1.1, 0.2),
            ylabel=r"$\bf{\%}$ $\bf{reads}$ $\bf{methylated}$",
            yaxis_formatter=mtick.PercentFormatter(1.0),
            xtick_rotation=0,
        )[0].save_organized(
            OUTPUT_DIR,
            "04-EMseq-controls",
            base + "-" + condition,
        )

# %% Create main Manhattan plots

importlib.reload(lib)

score_yticks = np.arange(-250, 251, 50)
score_yticklabels = abs(score_yticks)

for cell_line, _, _, base in comparisons_iter():
    targeted_gene = targeted_genes[cell_line]
    tss_index = lib.gene_info(
        gene_metadata,
        etest[base],
        targeted_gene,
    )["tss_index"]

    print(f"Working on {base} Manhattan plot...")

    lib.manhattan(
        etest[base].filter(pl.col("chr_order") < 25),
        by="chr",
        feature="score",
        two_sided=("targeting", "non-targeting"),
        highlight=None,
        yticks=score_yticks,
        yticklabels=score_yticklabels,
        ylabel=r"$\bf{Significance}$ $\bf{score}$" + "\n" + r"$-\log_{10}($FDR$)$",
        use_xticks=True,
        arrow=(targeted_gene, tss_index),
    )[0].save_organized(
        OUTPUT_DIR,
        "05-EMseq-manhattan",
        base,
        svg=False,
    )

# %% Create zoomed Manhattan plots

importlib.reload(lib)

effect_size_yticks = np.arange(-1, 1.1, 0.25)
effect_size_yticklabels = [round(yt, 2) for yt in effect_size_yticks]

for cell_line, _, _, base in comparisons_iter():
    targeted_gene = targeted_genes[cell_line]
    gi = lib.gene_info(gene_metadata, etest[base], targeted_gene)
    tss_index = gi["tss_index"]
    surround = gi["surround"]

    print(f"Working on zoomed {base} Manhattan plot...")

    lib.manhattan(
        etest[base].filter(surround),
        by="chr",
        feature="effect_size",
        two_sided=("targeting", "non-targeting"),
        highlight=None,
        yticks=effect_size_yticks,
        yticklabels=effect_size_yticklabels,
        ylabel=r"$\bf{Effect}$ $\bf{size}$",
        use_xticks=False,
        arrow=(targeted_gene, tss_index),
    )[0].save_organized(
        OUTPUT_DIR,
        "06-EMseq-zoomed-manhattan",
        base,
    )

# %% Create correlation plots

importlib.reload(lib)

for cell_line, _, _, base in comparisons_iter():
    data = (
        pl.read_csv(
            os.path.join(RNA_TEST_DIR, f"{base}.csv"),
            null_values="NA",
        )
        .filter(~pl.col("padj").is_null())
        .select(
            pl.col("").alias("ensembl_gene_id_version"),
            (
                pl.col("log2FoldChange").sign()
                * pl.col("padj")
                .log(
                    base=10,
                )
                .neg()
            ).alias("rscore"),
            pl.col("log2FoldChange").alias("expression_l2fc"),
            pl.col("padj").alias("expression_padj"),
        )
        .join(
            pl.read_csv(
                os.path.join(EM_AGG_DIR, f"{base}.csv"),
            ).select(
                pl.col("ensembl_gene_id_version"),
                (
                    pl.col("mean_diff").sign()
                    * pl.col("min_fdr")
                    .log(
                        base=10,
                    )
                    .neg()
                ).alias("escore"),
                pl.col("mean_diff").alias("methylation_diff"),
                pl.col("min_fdr").alias("methylation_fdr"),
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
        + r"$-\log_{10}($adjusted $p$-value$) \cdot $sign$($difference$)$",
        ylabel=r"$\bf{EMseq}$ $\bf{significance}$ $\bf{score}$"
        + "\n"
        + r"$-\log_{10}($FDR$) \cdot $sign$($difference$)$",
        highlight=pl.col("external_gene_name") == targeted_genes[cell_line],
    )[0].save_organized(
        OUTPUT_DIR,
        "07-correlation",
        base,
    )

    csv_dir = os.path.join(OUTPUT_DIR, "07-correlation", "csv")
    os.makedirs(csv_dir, exist_ok=True)

    export_data = data.select(
        "ensembl_gene_id_version",
        "external_gene_name",
        "expression_l2fc",
        "expression_padj",
        "methylation_diff",
        "methylation_fdr",
    )

    export_data.sort(by="expression_padj").head(10).write_csv(
        os.path.join(csv_dir, base + "-by-expression.csv"),
    )

    export_data.sort(by="methylation_fdr").head(10).write_csv(
        os.path.join(csv_dir, base + "-by-methylation.csv"),
    )

# %% Create methylation plots

importlib.reload(lib)

for cell_line, _, _, base in comparisons_iter():
    targeted_gene = targeted_genes[cell_line]
    gi = lib.gene_info(gene_metadata, etest[base], targeted_gene)
    tss = gi["tss"]
    surround = gi["surround"]

    lib.methylation_plot(
        etest[base].filter(surround),
        tracks=("non-targeting", "targeting"),
        targeted_gene=targeted_gene,
        tss=tss,
        surround=500,
    )[0].save_organized(
        OUTPUT_DIR,
        "08-EMseq-methylation",
        base,
    )
