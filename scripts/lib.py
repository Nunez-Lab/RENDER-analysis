import os

import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import polars as pl


def save(self, filename, exts=None, *args, **kwargs):
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)

    os.makedirs(dirname, exist_ok=True)

    if exts is not None:
        for ext in exts:
            new_dirname = os.path.join(dirname, ext)
            os.makedirs(new_dirname, exist_ok=True)
            self.savefig(
                os.path.join(new_dirname, basename + "." + ext),
                *args,
                **kwargs,
            )
    else:
        self.savefig(filename, *args, **kwargs)

    plt.close(self)


matplotlib.figure.Figure.save = save


def save_organized(self, output_dir, tag, base, sep="-", svg=True):
    path = os.path.join(output_dir, tag, tag + sep + base)
    if svg:
        self.save(path, exts=["png", "svg"])
    else:
        self.save(path, exts=["png"], dpi=300)


matplotlib.figure.Figure.save_organized = save_organized


def rna_count_plot(
    df,
    x_feature,
    y_feature,
    *,
    title,
    highlight=pl.lit(False),
    gene_name_feature="external_gene_name",
    xlabel=None,
    ylabel=None,
):
    xlabel = x_feature if xlabel is None else xlabel
    ylabel = y_feature if ylabel is None else ylabel

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    df = df.with_columns(
        (1 + df[x_feature]).log(base=2),
        (1 + df[y_feature]).log(base=2),
    )

    ax.scatter(
        df.filter(~highlight)[x_feature],
        df.filter(~highlight)[y_feature],
        c="0.5",
        zorder=5,
        alpha=0.5,
        s=5,
    )

    ax.scatter(
        df.filter(highlight)[x_feature],
        df.filter(highlight)[y_feature],
        c=PURPLE,
        zorder=10,
        alpha=1,
        s=10,
    )

    if gene_name_feature is not None:
        for row in df.filter(highlight).iter_rows(named=True):
            ax.annotate(
                text=row[gene_name_feature],
                xy=(row[x_feature], row[y_feature]),
                zorder=15,
                xytext=(3, 3),
                textcoords="offset pixels",
                ha="left",
                va="top",
                color=PURPLE,
            )

    ax.set_xlabel(
        r"$\bf{" + xlabel.replace("_", r"\_") + "}$\n" + r"$\log_{2}(1 + $TPM$)$"
    )

    ax.set_ylabel(
        r"$\bf{" + ylabel.replace("_", r"\_") + "}$\n" + r"$\log_{2}(1 + $TPM$)$"
    )

    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)

    ax.set_xticks(np.arange(0, 16.1, 2))
    ax.set_yticks(np.arange(0, 16.1, 2))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)

    fig.tight_layout()

    return fig, ax


def volcano_plot(
    df,
    *,
    title,
    treatment_name,
    control_name,
    highlight=pl.lit(False),
    threshold=None,
    show_threshold_line=False,
    gene_name_feature=None,
):
    if gene_name_feature is not None:
        assert threshold is not None

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    x = "log2FoldChange"
    y = "score"

    df = df.with_columns(
        **{
            y: -pl.col("padj").log(base=10),
        }
    ).filter(
        pl.col("baseMean") > 10,
    )

    ax.scatter(
        df.filter(~highlight)[x],
        df.filter(~highlight)[y],
        c="0.5",
        zorder=5,
        alpha=0.5,
        s=5,
    )

    ax.scatter(
        df.filter(highlight)[x],
        df.filter(highlight)[y],
        c=PURPLE,
        zorder=10,
        alpha=1,
        s=10,
    )

    if threshold is not None:
        if show_threshold_line:
            ax.axhline(
                y=threshold,
                ls="--",
                c="0.7",
                lw=1,
            )

        if gene_name_feature is not None:
            for row in df.filter(
                pl.col(y) > threshold,
                highlight,
            ).iter_rows(named=True):
                ax.annotate(
                    text=row[gene_name_feature],
                    xy=(row[x], row[y]),
                    zorder=15,
                    xytext=(3, 3),
                    textcoords="offset pixels",
                    ha="left",
                    va="top",
                    color=PURPLE,
                )

            for row in df.filter(
                pl.col(y) > threshold,
                ~highlight,
            ).iter_rows(named=True):
                ax.annotate(
                    text=row[gene_name_feature],
                    xy=(row[x], row[y]),
                    zorder=15,
                    xytext=(3, 3),
                    textcoords="offset pixels",
                    ha="left",
                    va="top",
                    color="0.5",
                )

    ax.set_xlabel(r"$\log_{2}($fold change$)$")
    ax.set_ylabel(r"$-\log_{10}($adjusted $p$-value$)$")

    xmax = df[x].abs().max() + 1
    ax.set_xlim(-xmax, xmax)

    ymax = df[y].max()
    if ymax < 30:
        ax.set_ylim(0, 30)
        ax.set_yticks(np.arange(0, 30.1, 5))
    else:
        ymax = max(150, ymax)
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(0, ymax + 0.1, 30))

    ax.text(
        0.03,
        1,
        f"Lower in\n{treatment_name}",
        color=PURPLE,
        fontweight="bold",
        ha="left",
        va="top",
        transform=ax.transAxes,
    )

    ax.text(
        1,
        1,
        f"Lower in\n{control_name}",
        color=BLUE,
        fontweight="bold",
        ha="right",
        va="top",
        transform=ax.transAxes,
    )

    ax.axvline(x=0, ls="--", lw=1, color="0.7")

    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    return fig, ax


def load_dss_results(path):
    return (
        pl.read_csv(path)
        .with_columns(
            mu_control=pl.col("mu1"),
            mu_treatment=pl.col("mu2"),
            effect_size=pl.col("diff").neg(),  # TODO swap order in DSS!
            score=pl.when(pl.col("mu1") > pl.col("mu2"))
            .then(pl.col("fdr").log10())
            .otherwise(-pl.col("fdr").log10()),
            chr=pl.when(pl.col("chr") == "phage_lambda")
            .then(pl.col("chr"))
            .when(pl.col("chr") == "plasmid_puc19c")
            .then(pl.col("chr"))
            .when(pl.col("chr").str.contains("_"))
            .then(pl.lit("other"))
            .otherwise(pl.col("chr")),
        )
        .with_columns(
            chr_order=pl.when(pl.col("chr") == "chrX")
            .then(23)
            .when(pl.col("chr") == "chrY")
            .then(24)
            .when(pl.col("chr") == "chrM")
            .then(25)
            .when(pl.col("chr") == "other")
            .then(26)
            .when(pl.col("chr") == "phage_lambda")
            .then(27)
            .when(pl.col("chr") == "plasmid_puc19c")
            .then(28)
            .otherwise(pl.col("chr").str.slice(3).str.to_integer(strict=False)),
            nice_chr=pl.when(pl.col("chr") == "phage_lambda")
            .then(pl.lit("Unmethylated control\n(lambda phage DNA)"))
            .when(pl.col("chr") == "plasmid_puc19c")
            .then(pl.lit("Methylated control\n(pUC19 plasmid DNA)"))
            .otherwise(pl.col("chr")),
            short_chr=pl.when(pl.col("chr") == "other")
            .then(pl.lit("?"))
            .when(pl.col("chr") == "phage_lambda")
            .then(pl.lit("Lambda"))
            .when(pl.col("chr") == "plasmid_puc19c")
            .then(pl.lit("pUC19"))
            .otherwise(pl.col("chr").str.slice(3)),
        )
        .select(
            "chr",
            "nice_chr",
            "short_chr",
            "chr_order",
            "pos",
            "mu_control",
            "mu_treatment",
            "score",
            "effect_size",
        )
        .sort(by=["chr_order", "pos"])
        .with_row_index()
    )


# Source: https://sronpersonalpages.nl/~pault/
BLUE = "#4477AA"
PURPLE = "#AA3377"
GREEN = "#228833"
RED = "#EE6677"


def manhattan(
    df,
    *,
    by,
    feature,
    two_sided,
    highlight,
    index="index",
    main_color=BLUE,
    secondary_color=PURPLE,
    highlight_color=GREEN,
    yticks=None,
    yticklabels=None,
    ylabel=None,
    yaxis_formatter=None,
    xtick_rotation=90,
    xlabel=None,
    figsize=(10, 3),
    use_xticks=True,
    arrow=None,
):
    if yticklabels is not None:
        assert yticks is not None

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    xticks = []
    xticklabels = []

    for (c,), g in df.group_by(by):
        lo = g[index].min()
        hi = g[index].max()
        mid = (lo + hi) / 2

        if use_xticks:
            xticks.append(mid)
            xticklabels.append(c)

        if two_sided is not None:
            g = g.with_columns(
                color=pl.when(pl.col(feature) > 0)
                .then(pl.lit(secondary_color))
                .otherwise(pl.lit(main_color))
            )
        else:
            g = g.with_columns(color=pl.lit(main_color))

        if highlight is not None:
            g = g.with_columns(
                color=pl.when(highlight)
                .then(pl.lit(highlight_color))
                .otherwise(pl.col("color"))
            )

        ax.scatter(
            g[index],
            g[feature].clip(upper_bound=yticks.max()),
            marker=",",
            s=1,
            c=g["color"],
        )

        if use_xticks:
            ax.axvline(
                lo,
                color="lightgray",
                linewidth=0.5,
            )

    left = df[index].min()
    right = df[index].max()

    if use_xticks:
        ax.axvline(right, color="lightgray", linewidth=0.5)

    ax.set_xlim(left, right)

    if xlabel:
        ax.set_xlabel(xlabel)

    ax.set_xticks(
        xticks,
        labels=xticklabels,
        rotation=xtick_rotation,
    )

    ax.set_ylabel(ylabel)

    if two_sided is not None:
        xpadding = 0.005
        ypadding = 0.02
        ax.text(
            1 - xpadding,
            1 - ypadding,
            f"Higher in {two_sided[0]}",
            color=secondary_color,
            fontweight="bold",
            ha="right",
            va="top",
            transform=ax.transAxes,
        )
        ax.text(
            1 - xpadding,
            0 + ypadding,
            f"Higher in {two_sided[1]}",
            color=main_color,
            fontweight="bold",
            ha="right",
            va="bottom",
            transform=ax.transAxes,
        )

    if arrow is not None:
        text = arrow[0]
        xy = (arrow[1], 1.01)
        xytext = (arrow[1], 1.12)

        ax.annotate(
            text,
            xy=xy,
            xytext=xytext,
            ha="center",
            va="bottom",
            fontweight="bold",
            xycoords=ax.get_xaxis_transform(),
            textcoords=ax.get_xaxis_transform(),
            arrowprops=dict(
                facecolor="black",
                headwidth=8,
                headlength=8,
            ),
        )

    if yticks is not None:
        ax.set_yticks(
            yticks,
            labels=yticklabels if yticklabels is not None else yticks,
        )
    elif two_sided is not None:
        ylim_lo, ylim_hi = ax.get_ylim()
        bound = max(abs(ylim_lo), abs(ylim_hi))
        ax.set_ylim(-bound, bound)

    if yaxis_formatter is not None:
        ax.yaxis.set_major_formatter(yaxis_formatter)

    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    return fig, ax


def gene_info(metadata, data, symbol):
    row = metadata.row(
        by_predicate=pl.col("external_gene_name") == symbol,
        named=True,
    )

    chr = "chr" + row["chromosome_name"]

    expr = (pl.col("chr") == chr) & (
        pl.col("pos").is_between(
            row["start_position"],
            row["end_position"],
        )
    )

    return {
        "tss_index": data.filter(expr)["index"].min(),
        "surround": (
            (pl.col("chr") == chr)
            & (
                pl.col("pos").is_between(
                    row["start_position"] - 10_000,
                    row["end_position"] + 10_000,
                )
            )
        ),
    }


def correlation_plot(
    df,
    s1,
    s2,
    *,
    xlabel,
    ylabel,
    gene_name_feature="external_gene_name",
    highlight=pl.lit(False),
):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    df = df.with_columns(
        **{
            s2: df[s2].clip(
                lower_bound=-250,
                upper_bound=250,
            )
        }
    )

    ax.scatter(
        df.filter(~highlight)[s1],
        df.filter(~highlight)[s2],
        c="0.5",
        zorder=5,
        alpha=0.5,
        s=5,
    )

    ax.scatter(
        df.filter(highlight)[s1],
        df.filter(highlight)[s2],
        c=PURPLE,
        zorder=10,
        alpha=1,
        s=10,
    )

    for row in df.filter(highlight).iter_rows(named=True):
        ax.annotate(
            text=row[gene_name_feature],
            xy=(row[s1], row[s2]),
            xytext=(3, 3),
            textcoords="offset pixels",
            ha="left",
            va="top",
            color=PURPLE,
            zorder=15,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    xmax = df.filter(pl.col(s1).is_finite())[s1].abs().max()
    if xmax < 30:
        xmax = 30
        ax.set_xlim(-xmax, xmax)
        ax.set_xticks(np.arange(-xmax, xmax + 0.1, 10))
    else:
        xmax = max(150, xmax)
        ax.set_xlim(-xmax, xmax)
        ax.set_yticks(np.arange(-xmax, xmax + 0.1, 30))

    ymax = max(250, df[s2].abs().max())
    ax.set_ylim(-ymax, ymax)
    ax.set_yticks(np.arange(-ymax, ymax + 0.1, 50))

    ax.axvline(x=0, lw=1, ls="--", color="0.8")
    ax.axhline(y=0, lw=1, ls="--", color="0.8")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    return fig, ax
