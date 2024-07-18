import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

#### Begining of scripts

import pandas as pd


def main(samples, completeness_files, bin_table):
    sample_data = {}
    div = {}

    df_list = []

    for i, sample in enumerate(samples):
        sample_data = pd.read_table(completeness_files[i], index_col=0)
        sample_data["Sample"] = sample

        df_list.append(sample_data)

    df = pd.concat(df_list, axis=0)

    df.to_csv(bin_table, sep="\t")


if __name__ == "__main__":
    main(
        samples=snakemake.params.samples,
        completeness_files=snakemake.input.completeness_files,
        bin_table=snakemake.output.bin_table,
    )
