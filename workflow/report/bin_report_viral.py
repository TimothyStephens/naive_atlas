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


from common_report import *

import pandas as pd
import plotly.express as px


from utils.taxonomy import tax2table


def make_plots(bin_info):
    div = {}

    div["input_file"] = f"{bin_info} and {snakemake.input.bins2species}"

    # Prepare data
    df = pd.read_table(bin_info, index_col=0)
    df["Bin Id"] = df.index  # need it also as column

    # add species info
    bin2species = pd.read_table(snakemake.input.bins2species, index_col=0)
    df = df.join(bin2species)

    logging.info(df.head())

    logging.info(bin2species.head())

    # calculate number of genomes/bins
    st = pd.DataFrame(columns=["Bins", "Species"])

    def add_stats(name, d):
        st.loc[name, "Bins"] = d.shape[0]
        st.loc[name, "Species"] = d.Representative.unique().shape[0]

    add_stats("All", df)

    df.eval("Quality_score = `virus_score`", inplace=True)
    div[
        "QualityScore"
    ] = "<p>Quality score is calculated as: virus_score.</p>"
    add_stats("Quality score >0.5 ", df.query("Quality_score>0.5"))
    add_stats("Good quality", df.query("`virus_score`>0.9 & `agreement`>0.5"))
    add_stats("Quality score >0.9 ", df.query("Quality_score>0.9"))

    div["table"] = st.to_html()

    logging.info(df.describe())

    # file             length  topology             coordinates  n_genes  genetic_code  virus_score  fdr     n_hallmarks  marker_enrichment  n_genes_with_taxonomy  agreement  taxid    lineage
    # format  type  num_seqs  sum_len  min_len  avg_len   max_len  Q1       Q2        Q3        sum_gap  N50      Q20(%)  Q30(%)  GC(%)
    hover_data = [
        "virus_score",
        "fdr",
        "n_hallmarks",
        "marker_enrichment",
        "n_genes_with_taxonomy",
        "agreement",
        "lineage",
        "N50",
        "GC(%)",
    ]
    size_name = "sum_len"

    lineage_name = "Species"

    # 2D plot

    logging.info("make 2d plot")
    fig = px.scatter(
        data_frame=df,
        y="virus_score",
        x="n_hallmarks",
        color=lineage_name,
        size=size_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(-0.1, 1.1))
    #fig.update_xaxes(range=(-0.2, 10.1))
    div["2D"] = fig.to_html(**HTML_PARAMS)

    # 2D plot

    logging.info("make 2d plot species")
    fig = px.scatter(
        data_frame=df.loc[df.Representative.unique()],
        y="virus_score",
        x="n_hallmarks",
        color=lineage_name,
        size=size_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(-0.1, 1.1))
    #fig.update_xaxes(range=(-0.2, 10.1))
    div["2Dsp"] = fig.to_html(**HTML_PARAMS)

    ## By sample
    logging.info("plot by sample")
    fig = px.strip(
        data_frame=df,
        y="Quality_score",
        x="Sample",
        color=lineage_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(-0.1, 1.1))
    div["bySample"] = fig.to_html(**HTML_PARAMS)

    # # By species
    # logging.info("plot by species")
    # fig = px.strip(
    #     data_frame=df,
    #     y="Quality_score",
    #     x=lineage_name,
    #     hover_data=hover_data,
    #     hover_name="Bin Id",
    # )
    # fig.update_yaxes(range=(50, 102))
    # div["byPhylum"] = fig.to_html(**HTML_PARAMS)

    return div


# main


div = make_plots(bin_info=snakemake.input.bin_info)


make_html(
    div=div,
    report_out=snakemake.output.report,
    html_template_file=os.path.join(reports_dir, "template_bin_report.html"),
    wildcards=snakemake.wildcards,
)
