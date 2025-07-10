################################################################################
# %% Imports

import importlib
import os

import matplotlib.ticker as mtick
import numpy as np
import polars as pl

try:
    import lib
except ModuleNotFoundError:
    import scripts.lib as lib

importlib.reload(lib)

QUICK = True

################################################################################
# %% Command-line arguments

# METADATA_DIR = sys.argv[1]
# GENE_METADATA_DIR = sys.argv[2]
#
# RNA_AGG_DIR = sys.argv[3]
# RNA_TEST_DIR = sys.argv[4]
#
# EM_AVG_PATH = sys.argv[5]
# EM_TEST_DIR = sys.argv[6]
# EM_AGG_DIR = sys.argv[7]
#
# OUTPUT_DIR = sys.argv[8]

################################################################################
# %% RNA-seq count plots

METADATA_DIR = "metadata"
GENE_METADATA_DIR = "output/gene-metadata"

RNA_AGG_DIR = "output/RNAseq/aggregated-reads"
RNA_TEST_DIR = "output/RNAseq/deseq2"

EM_AVG_PATH = "output/EMseq/average-methylation-info/info.tsv"
EM_TEST_DIR = "output/EMseq/dss"
EM_AGG_DIR = "output/EMseq/dss-aggregated"

OUTPUT_DIR = "output/analysis"

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
    pl.read_csv(os.path.join(RNA_AGG_DIR, "abundance.csv"))
    .rename({"": "ensembl_gene_id_version"})
    .join(
        gene_metadata,
        on="ensembl_gene_id_version",
        how="left",
        validate="1:1",
    )
)

rcounts = (
    pl.read_csv(os.path.join(RNA_AGG_DIR, "counts.csv"))
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
    if control == "untreated":
        continue
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
        highlight=pl.col("external_gene_name") == targeted_genes[cell_line],
        xlabel="Replicate\\ 1",
        ylabel="Replicate\\ 2",
        show_r2=True,
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

    if control == "untreated":
        xlabel = "Untreated"
    elif control == "non_targeting":
        xlabel = r"Non\!-\!targeting"
    else:
        raise ValueError(f"unknown control '{control}'")

    lib.rna_count_plot(
        df,
        "control",
        "treatment",
        xlabel=xlabel,
        ylabel="Treated",
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
        .with_columns(score=-pl.col("padj").log(base=10))
    )

    csv_dir = os.path.join(OUTPUT_DIR, "03-RNAseq-volcano", "csv")
    os.makedirs(csv_dir, exist_ok=True)

    rtest.filter(~pl.col("padj").is_null()).sort(by="padj").write_csv(
        os.path.join(csv_dir, base + ".csv"),
    )

    for kind, df in [("ABUNDANCE", rabundance), ("COUNTS", rcounts)]:
        rtest.filter(
            pl.col("score").abs() > lib.SCORE_THRESHOLD,
            pl.col("log2FoldChange").abs() > lib.L2FC_THRESHOLD,
        ).join(
            df,
            on="ensembl_gene_id_version",
            how="left",
            validate="1:1",
        ).sort(by="padj").write_csv(
            os.path.join(csv_dir, f"HITS-{kind}-{base}.csv"),
        )

    if control == "untreated":
        control_name = "untreated"
    elif control == "non_targeting":
        control_name = "non-targeting"
    else:
        raise ValueError(f"unknown control '{control}'")

    lib.volcano_plot(
        rtest,
        title=f"RENDER {treatment} vs {control} for {cell_line}",
        treatment_name="treated",
        control_name=control_name,
        highlight=pl.col("external_gene_name") == targeted_genes[cell_line],
        gene_name_feature="external_gene_name",
        threshold=False,
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

if not QUICK:
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

for cell_line, control, _, base in comparisons_iter():
    if control == "untreated":
        continue

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
        xlabel=r"$\bf{Differential}$ $\bf{expression}$"
        + "\n"
        + r"$-\log_{10}($adjusted $p$-value$) \cdot $sign$($difference$)$",
        ylabel=r"$\bf{Differential}$ $\bf{methylation}$"
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

    export_data.sort(by="expression_padj").write_csv(
        os.path.join(csv_dir, base + "-by-expression.csv"),
    )

    export_data.sort(by="methylation_fdr").write_csv(
        os.path.join(csv_dir, base + "-by-methylation.csv"),
    )

# %% Create methylation plots

importlib.reload(lib)

for cell_line, _, _, base in comparisons_iter():
    if control == "untreated":
        continue

    targeted_gene = targeted_genes[cell_line]
    surround_amount = 2000
    gi = lib.gene_info(
        gene_metadata,
        etest[base],
        targeted_gene,
        surround_amount=surround_amount,
    )
    tss = gi["tss"]
    surround = gi["surround"]

    lib.methylation_plot(
        etest[base].filter(surround),
        tracks=("Non-targeting", "Targeting"),
        targeted_gene=targeted_gene,
        tss=tss,
        surround=surround_amount,
    )[0].save_organized(
        OUTPUT_DIR,
        "08-EMseq-methylation",
        base,
    )
