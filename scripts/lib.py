import os

import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import polars as pl
import scipy.stats

# Source: https://sronpersonalpages.nl/~pault/
BLUE = "#4477AA"
PURPLE = "#AA3377"
GREEN = "#228833"
PINK = "#EE6677"
RED = "#880000"
YELLOW = "#CCBB44"

NON_HIGHLIGHTED_COLOR = "0.6"
HIGHLIGHT_FONT_OFFSET = 4

L2FC_THRESHOLD = 2
SCORE_THRESHOLD = 5

plt.rcParams["font.family"] = "Arial"


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
    annotation_config=("left", "bottom", 3, 0),  # (ha, va, x offset, y offset)
    title=None,
    highlight=pl.lit(False),
    gene_name_feature="external_gene_name",
    xlabel=None,
    ylabel=None,
    highlight_color=RED,
    fontsize=13,
    dotsize=25,
    show_r2=False,
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
        c=NON_HIGHLIGHTED_COLOR,
        zorder=5,
        alpha=0.5,
        s=dotsize,
    )

    ax.scatter(
        df.filter(highlight)[x_feature],
        df.filter(highlight)[y_feature],
        c=highlight_color,
        zorder=10,
        alpha=1,
        s=dotsize,
    )

    if gene_name_feature is not None:
        for row in df.filter(highlight).iter_rows(named=True):
            ax.annotate(
                text=row[gene_name_feature],
                xy=(row[x_feature], row[y_feature]),
                zorder=15,
                xytext=(annotation_config[2], annotation_config[3]),
                textcoords="offset pixels",
                ha=annotation_config[0],
                va=annotation_config[1],
                color=highlight_color,
                fontsize=fontsize + HIGHLIGHT_FONT_OFFSET,
            )

    if show_r2:
        result = scipy.stats.pearsonr(df[x_feature], df[y_feature])
        r2 = result.statistic**2
        ax.text(
            0.05,
            0.95,
            "$R^2 = " + str(round(r2, 3)) + "$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=fontsize,
        )

    ax.set_xlabel(
        r"$\bf{" + xlabel.replace("_", r"\_") + "}$\n" + r"$\log_{2}(1 + $TPM$)$"
    )

    ax.set_ylabel(
        r"$\bf{" + ylabel.replace("_", r"\_") + "}$\n" + r"$\log_{2}(1 + $TPM$)$"
    )

    if title is not None:
        ax.set_title(title)

    ax.spines[["top", "right"]].set_visible(False)

    lim = 16  # max(df[x_feature].max(), df[y_feature].max(), 16)
    ax.set_xticks(np.arange(0, lim + 0.1, 2))
    ax.set_yticks(np.arange(0, lim + 0.1, 2))
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)

    fig.tight_layout()

    return fig, ax


def rscore_range(maxval):
    extra = 50 if int(maxval) % 100 > 90 else 0
    return np.ceil(maxval / 50) * 50 + extra


def volcano_plot(
    df,
    *,
    title,
    treatment_name,
    control_name,
    highlight=pl.lit(False),
    gene_name_feature=None,
    highlight_color=RED,
    fontsize=13,
    dotsize=25,
    threshold=True,
):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    x = "log2FoldChange"
    y = "score"

    ax.scatter(
        df.filter(~highlight)[x],
        df.filter(~highlight)[y],
        c=NON_HIGHLIGHTED_COLOR,
        zorder=5,
        alpha=0.5,
        s=dotsize,
    )

    ax.scatter(
        df.filter(highlight)[x],
        df.filter(highlight)[y],
        c=highlight_color,
        zorder=10,
        alpha=1,
        s=dotsize,
    )

    if gene_name_feature is not None:
        for row in df.filter(highlight).iter_rows(named=True):
            ax.annotate(
                text=row[gene_name_feature],
                xy=(row[x], row[y]),
                zorder=15,
                xytext=(-3, -3),
                textcoords="offset pixels",
                ha="right",
                va="top",
                color=highlight_color,
                fontsize=fontsize + HIGHLIGHT_FONT_OFFSET,
            )

    ax.set_xlabel(r"$\log_{2}($fold change$)$")
    ax.set_ylabel(r"$-\log_{10}($adjusted $p$-value$)$")

    xmax = max(df[x].abs().max(), 16)
    ax.set_xlim(-xmax, xmax)
    ax.set_xticks(np.arange(-xmax, xmax + 0.1, 4))

    ymax = rscore_range(df[y].max())
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax + 0.1, 50))

    ax.text(
        0.03,
        1,
        f"Lower in\n{treatment_name}",
        color=PURPLE,
        fontweight="bold",
        ha="left",
        va="top",
        transform=ax.transAxes,
        fontsize=fontsize,
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
        fontsize=fontsize,
    )

    if threshold:
        ax.axvline(x=-L2FC_THRESHOLD, ls="--", lw=1, color="0.8")
        ax.axvline(x=+L2FC_THRESHOLD, ls="--", lw=1, color="0.8")
        ax.axhline(y=SCORE_THRESHOLD, ls="--", lw=1, color="0.8")

    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    return fig, ax


def load_dss_results(path):
    return (
        pl.read_csv(path)
        .with_columns(
            mu_treatment=pl.col("mu1"),
            mu_control=pl.col("mu2"),
            effect_size=pl.col("diff"),
            score=pl.col("fdr").log10().neg() * pl.col("diff").sign(),
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


def gene_info(metadata, data, symbol, surround_amount=10_000):
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
        "tss": row["start_position"],
        "tss_index": data.filter(expr)["index"].min(),
        "surround": (
            (pl.col("chr") == chr)
            & (
                pl.col("pos").is_between(
                    row["start_position"] - surround_amount,
                    row["end_position"] + surround_amount,
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
    highlight_color=RED,
    fontsize=13,
    dotsize=25,
    s2_bound=200,
):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    df = df.with_columns(
        **{
            s2: df[s2].clip(
                lower_bound=-s2_bound,
                upper_bound=s2_bound,
            )
        }
    )

    ax.scatter(
        df.filter(~highlight)[s1],
        df.filter(~highlight)[s2],
        c=NON_HIGHLIGHTED_COLOR,
        zorder=5,
        alpha=0.5,
        s=dotsize,
    )

    ax.scatter(
        df.filter(highlight)[s1],
        df.filter(highlight)[s2],
        c=highlight_color,
        zorder=10,
        alpha=1,
        s=dotsize,
    )

    for row in df.filter(highlight).iter_rows(named=True):
        ax.annotate(
            text=row[gene_name_feature],
            xy=(row[s1], row[s2]),
            xytext=(3, -3),
            textcoords="offset pixels",
            ha="left",
            va="top",
            color=highlight_color,
            zorder=15,
            fontsize=fontsize + HIGHLIGHT_FONT_OFFSET,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    xmax = df[s1].abs().max()
    if xmax < 200:
        xmax = 200
    else:
        xmax = max(xmax, 400)
    ax.set_xlim(-xmax, xmax)
    ax.set_xticks(np.arange(-xmax, xmax + 0.1, 100))

    ymax = max(s2_bound, df[s2].abs().max())
    ax.set_ylim(-ymax, ymax)
    ax.set_yticks(np.arange(-ymax, ymax + 0.1, 50))

    ax.axvline(x=-SCORE_THRESHOLD, lw=1, ls="--", color="0.8")
    ax.axvline(x=+SCORE_THRESHOLD, lw=1, ls="--", color="0.8")

    ax.axhline(y=-SCORE_THRESHOLD, lw=1, ls="--", color="0.8")
    ax.axhline(y=+SCORE_THRESHOLD, lw=1, ls="--", color="0.8")

    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    return fig, ax


def methylation_plot(
    df,
    *,
    tracks,
    targeted_gene,
    tss,
    surround,
    figsize=(10, 3),
):
    fig, ax = plt.subplots(
        2,
        1,
        figsize=figsize,
        sharex=True,
    )

    conditions = ["control", "treatment"]
    colors = [BLUE, PURPLE]

    for i in range(len(conditions)):
        condition = conditions[i]
        color = colors[i]

        ax[i].vlines(
            df["pos"],
            0,
            1,
            colors="0.9",
            zorder=1,
            linewidth=1,
        )

        ax[i].vlines(
            df["pos"],
            0,
            df[f"mu_{condition}"],
            color=color,
            zorder=2,
            linewidth=1,
        )

        # ax[i].bar(
        #    df["pos"],
        #    np.ones_like(df["pos"]),
        #    width=3,
        #    color="0.9",
        #    zorder=1,
        #    hatch="----",
        #    edgecolor="white",
        #    linewidth=0,
        # )

        # ax[i].bar(
        #     df["pos"],
        #     df[f"mu_{condition}"],
        #     width=3,
        #     color=color,
        #     zorder=2,
        # )

        ax[i].spines[["top", "right"]].set_visible(False)

        ax[i].yaxis.set_label_position("right")

        ax[i].set_ylabel(
            tracks[i],
            color=colors[i],
            fontweight="bold",
            rotation=270,
            fontsize=11,
            labelpad=15,
        )

        ax[i].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

        ax[i].set_ylim(0, 1)

    lo = tss - surround
    hi = tss + surround

    ax[0].set_xlim(lo, hi)

    xticks = np.arange(lo, hi + 0.1, 500)
    offsets = xticks - tss
    labels = []
    for off in offsets:
        if off > 0:
            label = "+" + str(round(off)) + " bp"
        elif off == 0:
            label = ""
        else:
            label = str(round(off)) + " bp"
        labels.append(label)

    ax[0].set_xticks(xticks, labels=labels)

    ax[1].text(
        tss,
        -0.09,
        targeted_gene + " TSS",
        ha="center",
        va="top",
        fontweight="bold",
    )

    #  ax[0].annotate(
    #      targeted_gene,
    #      xy=(tss, 1.01),
    #      xytext=(tss, 1.27),
    #      ha="center",
    #      va="bottom",
    #      fontweight="bold",
    #      xycoords=ax[0].get_xaxis_transform(),
    #      textcoords=ax[0].get_xaxis_transform(),
    #      arrowprops=dict(
    #          facecolor="black",
    #          headwidth=8,
    #          headlength=8,
    #      ),
    #  )

    fig.supylabel(
        "% reads methylated",
        fontweight="bold",
        fontsize=13,
    )

    fig.tight_layout()

    return fig, ax
