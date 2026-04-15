#!/usr/bin/env python
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


# begining of script

import datetime
import shutil
import os
from snakemake.shell import shell


timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%X")


def calculate_read_length_hist(params_in, params_out):
    "get read lengths by running readlength.sh"
    if not params_in is None:
        shell(
            f" readlength.sh {params_in} out={params_out}"
            f" 1>>{snakemake.log[0]} 2>&1 "
        )
    else:
        shell(f"echo -e '#Reads:\t0\n#Bases:\t0\n#Max:\t0\n#Min:\t0\n#Avg:\t0.0\n#Median:\t0\n#Mode:\t0\n#Std_Dev:\t0.0\n' > {params_out}")


def calculate_insert_size_hist(params_in, params_out=None):
    "get read insert size by running bbmerge.sh"
    if not params_in is None:
        shell(
            f"bbmerge.sh "
            f" -Xmx{snakemake.resources.java_mem}M "
            f" threads={snakemake.threads} "
            f" {params_in} "
            f" {snakemake.params.flags} k={snakemake.params.kmer} "
            f" extend2={snakemake.params.extend2} "
            f" ihist={params_out} merge=f "
            f" mininsert0=35 minoverlap0=8 "
            f" prealloc=t prefilter=t "
            f" minprob={snakemake.params.minprob} "
            f" 1>>{snakemake.log[0]} 2>&1"
        )
    else:
        shell(f"echo -e '#Mean\t0.0\n#Median\t0\n#Mode\t0\n#STDev\t0.0\n#PercentOfPairs\t0.0\n' > {params_out}")


# Generate stats for each read type if present.
if hasattr(snakemake.input, 'R1'):
    calculate_read_length_hist(f"in1={snakemake.input.R1} in2={snakemake.input.R2}", f"{snakemake.output.lenHist_pe}")
    calculate_insert_size_hist(f"in1={snakemake.input.R1} in2={snakemake.input.R2}", f"{snakemake.output.insertHist_pe}")
else:
    calculate_read_length_hist(None,                                                 f"{snakemake.output.lenHist_pe}")
    calculate_insert_size_hist(None,                                                 f"{snakemake.output.insertHist_pe}")


if hasattr(snakemake.input, 'SE'):
    calculate_read_length_hist(f"in={snakemake.input.SE}", f"{snakemake.output.lenHist_se}")
    calculate_insert_size_hist(None,                       f"{snakemake.output.insertHist_se}")
else:
    calculate_read_length_hist(None,                       f"{snakemake.output.lenHist_se}")
    calculate_insert_size_hist(None,                       f"{snakemake.output.insertHist_se}")


if hasattr(snakemake.input, 'LR'):
    calculate_read_length_hist(f"in={snakemake.input.LR}", f"{snakemake.output.lenHist_lr}")
    calculate_insert_size_hist(None,                       f"{snakemake.output.insertHist_lr}")
else:
    calculate_read_length_hist(None,                       f"{snakemake.output.lenHist_lr}")
    calculate_insert_size_hist(None,                       f"{snakemake.output.insertHist_lr}")
