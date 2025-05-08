import polars as pl
import matplotlib.pyplot as plt

OUTPUT_DIR = "output/RNA_seq_Replicate_2"

def preprocess(df):
    return df.filter(pl.col("target_id").str.starts_with("NM_"))

def rpk(df):
    return df.with_columns(
        rpk=pl.col("est_counts") / (pl.col("eff_length") / 1000)
    )

def combine(df1, df2, suffix):
    df1 = rpk(df1)
    df2 = rpk(df2)
    return df1.join(
        df2,
        on="target_id",
        suffix="2",
        validate="1:1"
    ).with_columns(
        total_rpk=pl.col("rpk") + pl.col("rpk2"),
    ).with_columns(
        tpm1=pl.col("rpk") / (pl.col("rpk").sum() / 1_000_000),
        tpm2=pl.col("rpk2") / (pl.col("rpk2").sum() / 1_000_000),
        total_tpm=pl.col("total_rpk") / (pl.col("total_rpk").sum() / 1_000_000),
    ).select(
        pl.col("target_id"),
        pl.col("total_tpm").alias("tpm_" + suffix),
        pl.col("total_tpm").log(base=2).alias("log_tpm_" + suffix),
        pl.col("tpm1").alias("tpm_" + suffix + "1"),
        pl.col("tpm1").log(base=2).alias("log_tpm_" + suffix + "1"),
        pl.col("tpm2").alias("tpm_" + suffix + "2"),
        pl.col("tpm2").log(base=2).alias("log_tpm_" + suffix + "2"),
    )

# Load in the abundance data table (separated by tabs, not commas)
n1 = preprocess(pl.read_csv(f"{OUTPUT_DIR}/quant/DX5N1_S8/abundance.tsv", separator="\t"))
n2 = preprocess(pl.read_csv(f"{OUTPUT_DIR}/quant/DX5N2_S21/abundance.tsv", separator="\t"))
p1 = preprocess(pl.read_csv(f"{OUTPUT_DIR}/quant/DX5P1_S10/abundance.tsv", separator="\t"))
p2 = preprocess(pl.read_csv(f"{OUTPUT_DIR}/quant/DX5P2_S23/abundance.tsv", separator="\t"))

n = combine(n1, n2, suffix="n")
p = combine(p1, p2, suffix="p")

CD55_IDS = [
    "NM_000574.5",
    "NM_001114752.3",
    "NM_001300902.2",
    "NM_001300903.2",
    "NM_001300904.2",
]

data = n.join(p, on="target_id").with_columns(
    color=pl.when(pl.col("target_id").is_in(CD55_IDS)).then(pl.lit("red")).otherwise(pl.lit("gray"))
).filter(pl.col("tpm_n") > 0.05)

fig, ax = plt.subplots(1, 1, figsize=(4, 4))

for (c,), g in data.group_by("color"):
    if c == "red":
        alpha = 1
        zorder = 10
    else:
        alpha = 0.1
        zorder = 5
    ax.scatter(g["log_tpm_n1"], g["log_tpm_n2"], c=c, alpha=alpha, zorder=zorder)

ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)

fig.tight_layout()

fig.savefig(OUTPUT_DIR + "/plot.png")
