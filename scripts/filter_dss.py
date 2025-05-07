import sys
import os

import polars as pl


def usage():
    print(
        "usage: {sys.argv[0]} <minimum read count for each sample> <path to combined files> <desired path to filtered output>"
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    if len(sys.argv) != 4:
        print("error: incorrect number of arguments")
        usage()
        sys.exit(1)

    MIN_READ_COUNT = int(sys.argv[1])
    INPUT_PATH = sys.argv[2]
    OUTPUT_PATH = sys.argv[3]

    for filename in os.listdir(INPUT_PATH):
        pl.read_csv(
            f"{INPUT_PATH}/{filename}",
            separator="\t",
        ).filter(
            pl.col("N") >= MIN_READ_COUNT,
        ).write_csv(
            f"{OUTPUT_PATH}/{filename}",
            separator="\t",
        )
