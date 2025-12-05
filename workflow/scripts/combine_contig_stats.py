import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


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


import pandas as pd
from utils.parsers_samtools import parse_samtools_stats_file


def parse_map_stats(sample_data, out_tsv):
    sample_stats = {}
    for sample in sample_data.keys():
        df = pd.read_csv(sample_data[sample]["contig_stats"], sep="\t")

        assert df.shape[0] == 1, "Assumed only one row in file {}; found {}".format(
            sample_data[sample]["contig_stats"], df.iloc[0]
        )

        # n genes
        genes_df = pd.read_csv(sample_data[sample]["gene_table"], index_col=0, sep="\t")
        df["N_Predicted_Genes"] = genes_df.shape[0]

        # mappingt stats
        mapping_stats = parse_samtools_stats_file(sample_data[sample]["samtools_stats"])
        df["Total_Reads"]             = int(mapping_stats["raw total sequences"])
        df["Mapped_Reads"]            = int(mapping_stats["reads mapped"])
        df["Percent_Assembled_Reads"] = round((df["Mapped_Reads"]/df["Total_Reads"])*100, 2)

        logging.info(f"Stats for sample {sample}\n{df}")

        sample_stats[sample] = df

    stats_df = pd.concat(sample_stats, axis=0)
    stats_df.index = stats_df.index.get_level_values(0)
    # remove contig stats and keep only scaffold stats
    stats_df = stats_df.loc[:, ~stats_df.columns.str.startswith("scaf_")]
    stats_df.columns = stats_df.columns.str.replace("ctg_", "")
    # save
    stats_df.to_csv(out_tsv, sep="\t")
    return stats_df


def main(samples, contig_stats, gene_tables, samtools_stats, combined_stats):
    sample_data = {}
    for sample in samples:
        sample_data[sample] = {}
        for c_stat in contig_stats:
            # underscore version was for simplified local testing
            # if "%s_" % sample in c_stat:
            if "%s/" % sample in c_stat:
                sample_data[sample]["contig_stats"] = c_stat
        for g_table in gene_tables:
            # if "%s_" % sample in g_table:
            if "%s/" % sample in g_table:
                sample_data[sample]["gene_table"] = g_table
        for samtools_stat in samtools_stats:
            # if "%s_" % sample in samtools_stat:
            if "%s/" % sample in samtools_stat:
                sample_data[sample]["samtools_stats"] = samtools_stat

    parse_map_stats(sample_data, combined_stats)


if __name__ == "__main__":
    main(
        samples=snakemake.params.samples,
        contig_stats=snakemake.input.contig_stats,
        gene_tables=snakemake.input.gene_tables,
        samtools_stats=snakemake.input.samtools_stats,
        combined_stats=snakemake.output.combined_contig_stats,
    )
