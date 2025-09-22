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

#### Begining of scripts

from common_report import *


import pandas as pd
import plotly.express as px
from plotly import subplots
import plotly.graph_objs as go
import numpy as np


labels = {"Total_Reads": "Total Reads", "Total_Bases": "Total Bases"}


PLOT_PARAMS = dict(labels=labels)


import zipfile


def get_stats_from_zips(zips, samples):
    # def get_read_stats(samples, step):
    quality_pe = pd.DataFrame()
    quality_se = pd.DataFrame()
    quality_lr = pd.DataFrame()
    for zfile, sample in zip(zips, samples):
        zf = zipfile.ZipFile(zfile)

        if "pe/boxplot_quality.txt" in zf.namelist():
            with zf.open("pe/boxplot_quality.txt") as f:
                df = pd.read_csv(f, index_col=0, sep="\t")
                df.columns = [df.columns, [sample] * df.shape[1]]
                quality_pe = pd.concat(
                    (quality_pe, df[["mean_1", "mean_2"]]), axis=1
                )
        if "se/boxplot_quality.txt" in zf.namelist():
            with zf.open("se/boxplot_quality.txt") as f:
                df = pd.read_csv(f, index_col=0, sep="\t")
                quality_se[sample] = df.mean_1
        if "lr/boxplot_quality.txt" in zf.namelist():
            with zf.open("lr/boxplot_quality.txt") as f:
                df = pd.read_csv(f, index_col=0, sep="\t")
                quality_lr[sample] = df.mean_1

    return quality_pe, quality_se, quality_lr


def get_pe_read_quality_plot(df, quality_range, color_range):
    fig = subplots.make_subplots(cols=2)

    for i, sample in enumerate(df["mean_1"].columns):
        fig.append_trace(
            go.Scatter(
                x=df.index,
                y=df["mean_1"][sample].values,
                type="scatter",
                name=sample,
                legendgroup=sample,
                marker=dict(color=color_range[i]),
            ),
            1,
            1,
        )

        fig.append_trace(
            dict(
                x=df.index,
                y=df["mean_2"][sample].values,
                type="scatter",
                name=sample,
                legendgroup=sample,
                showlegend=False,
                marker=dict(color=color_range[i]),
            ),
            1,
            2,
        )

    fig.update_layout(
        yaxis=dict(range=quality_range, autorange=True, title="Average quality score"),
        xaxis1=dict(title="Position forward read"),
        xaxis2=dict(autorange="reversed", title="Position reverse read"),
    )

    return fig


def draw_se_read_quality(df, quality_range, color_range):
    fig = subplots.make_subplots(cols=1)

    for i, sample in enumerate(df.columns):
        fig.append_trace(
            go.Scatter(
                x=df.index,
                y=df[sample].values,
                type="scatter",
                name=sample,
                legendgroup=sample,
                marker=dict(color=color_range[i]),
            ),
            1,
            1,
        )

    fig.update_layout(
        yaxis=dict(range=quality_range, autorange=True, title="Average quality score"),
        xaxis=dict(title="Position read"),
    )
    return fig


def make_plots(
    samples, zipfiles_QC, read_counts, read_length, min_quality, insert_size_stats
):
    div = {}
    df = pd.read_csv(read_counts, index_col=[0, 1], sep="\t")
    
    
    ##
    ## Number of reads that went through the quality control process (per step)
    ##
    total_reads = df.Total_Reads.unstack()
    
    fig = px.bar(
        data_frame=total_reads,
        barmode="group",
        labels={"value": "Reads"},
        category_orders={"Step": ["1_raw", "2_deduplicated", "3_quality_filtered", "4_decontaminated", "5_final"]},
    )
    
    fig.update_yaxes(title="Number of reads")
    fig.update_xaxes(tickangle=45)
    # fig.update_layout(hovermode="x unified")
    
    div["Reads"] = fig.to_html(**HTML_PARAMS)
    
    
    
    ##
    ## Total (and per read type) number of reads/bases after QC
    ##
    data_qc = df.query('Step=="5_final"').reset_index()
    
    for var in ["Total_Reads", "Total_Bases", "Reads_pe", "Bases_pe", "Reads_se", "Bases_se", "Reads_lr", "Bases_lr"]:
        fig = px.strip(data_qc, y=var, hover_data=["Sample", var], **PLOT_PARAMS)
        fig.update_yaxes(range=(0, data_qc[var].max() * 1.1))
        div[var] = fig.to_html(**HTML_PARAMS)
    
    
    
    ##
    ## Quality values (per read type) along read
    ##
    N = len(samples)
    color_range = [
        "hsl(" + str(h) + ",50%" + ",50%)" for h in np.linspace(0, 360, N + 1)
    ]
    
    # load quality profiles for QC and low
    Quality_QC_pe, Quality_QC_se, Quality_QC_lr = get_stats_from_zips(zipfiles_QC, samples)
    
    # detrmine range of quality values and if paired
    max_quality = 1 + np.nanmax((Quality_QC_pe.max().max(), Quality_QC_se.max().max(), Quality_QC_lr.max().max()))
    quality_range = [min_quality, max_quality]
    
    # create plots if paired or not

    if Quality_QC_pe.shape[0] > 0:
        div["quality_QC_PE"] = get_pe_read_quality_plot(
            Quality_QC_pe, quality_range, color_range
        ).to_html(**HTML_PARAMS)
    else:
        div["quality_QC_PE"] = (
            "<p>No PE reads found in any of the samples.</p>"
        )

    if Quality_QC_se.shape[0] > 0:
        div["quality_QC_SE"] = draw_se_read_quality(
            Quality_QC_se, quality_range, color_range
        ).to_html(**HTML_PARAMS)
    else:
        div["quality_QC_SE"] = (
            "<p>No SE reads found in any of the samples.</p>"
        )

    if Quality_QC_lr.shape[0] > 0:
        # Use SE plot function since they are functionally the same
        div["quality_QC_LR"] = draw_se_read_quality(
            Quality_QC_lr, quality_range, color_range
        ).to_html(**HTML_PARAMS)
    else:
        div["quality_QC_LR"] = (
            "<p>No LR reads found in any of the samples.</p>"
        )
    
    
    
    ##
    ## Read length (per read type) plots
    ##
    data_length = pd.read_table(read_length, index_col=0).T
    data_length.index.name = "Sample"
    
    for fraction in ["PE", "SE", "LR"]:
        fig = px.bar(
            data_frame=data_length[data_length['fraction'].str.contains(fraction)],
            x="Median",
            error_x="Max",
            error_x_minus="Min",
            hover_data=["Median", "Max", "Min", "Avg", "Std_Dev", "Mode"],
        )
        
        fig.update_xaxes(title=fraction + " Read length")
        
        div["Length_"+fraction] = fig.to_html(**HTML_PARAMS)
    
    
    
    ##
    ## Insert size (PE only)
    ##
    data_insert = pd.read_table(insert_size_stats, index_col=0)
    data_insert.index.name = "Sample"
    
    fig = px.bar(
        data_frame=data_insert,
        x="Mean",
        error_x="STDev",
        hover_data=["Mean", "Median", "Mode", "PercentOfPairs"],
        labels={"PercentOfPairs": "Percent of pairs"},
    )
    
    fig.update_xaxes(title="Insert size")
    
    div["Insert"] = fig.to_html(**HTML_PARAMS)
    
    
    ## Return plots
    return div



# Make HTML report
div = make_plots(
    samples=snakemake.params.samples,
    zipfiles_QC=snakemake.input.zipfiles_QC,
    read_counts=snakemake.input.read_counts,
    read_length=snakemake.input.read_length_stats,
    min_quality=snakemake.params.min_quality,
    insert_size_stats=snakemake.input.read_insert_size_stats,
)

make_html(
    div=div,
    report_out=snakemake.output.report,
    html_template_file=os.path.join(reports_dir, "template_QC_report.html"),
)
