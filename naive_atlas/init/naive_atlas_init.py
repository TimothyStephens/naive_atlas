import os, sys
import logging
from ..color_logger import logger
import pandas as pd
import numpy as np
import click
from pathlib import Path
from .create_config import create_config
from .create_sample_table import get_samples_from_fastq, simplify_sample_names
from ..sample_table import (
    validate_sample_table,
    validate_bingroup_size,
    BinGroupSizeError,
)


from ..logo import print_logo
class LogoCommand(click.Command):
    """Click Command subclass that prints the logo before the help message."""
    def format_help(self, ctx, formatter):
        print_logo()
        super().format_help(ctx, formatter)



def prepare_sample_table_for_atlas(
    sample_table, assembler, outfile="samples.tsv"
):
    """
    Write the file `samples.tsv` and complete the sample names and paths for all
    files in `path`.
    Args:
            path_to_fastq (str): fastq/fasta data directory
    """

    if os.path.exists(outfile):
        logger.error(
            f"Output file {outfile} already exists I don't dare to overwrite it."
        )
        exit(1)

    simplify_sample_names(sample_table)

    sample_table["Bin_group"] = "all"
    sample_table["Assembler"] = assembler

    validate_sample_table(sample_table)

    sample_table.to_csv(outfile, sep="\t", index_label="ID")



@click.command(
    "init",
    cls=LogoCommand,
    short_help="Prepare configuration file and sample table for atlas run",
)
@click.argument("path_to_fastq", type=click.Path(readable=True))
@click.option(
    "-d",
    "--db-dir",
    default=os.path.join(os.path.realpath("."), "databases"),
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="Location to store databases (need ~150GB)",
)
@click.option(
    "-t",
    "--temp-dir",
    default=os.path.join(os.path.realpath("."), "tmp"),
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="Location to store temp files",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    default=".",
    help="Location to run naive_atlas",
)
@click.option(
    "-a",
    "--assembler",
    type=str,
    show_default=True,
    default="spades",
    help="Assembler to use",
)
@click.option(
    "--logger-debug/--logger-no-debug",
    default=False,
    show_default=True,
    help="Set logger level to debug.",
)
def run_init(
    path_to_fastq,
    db_dir,
    temp_dir,
    working_dir,
    assembler,
    logger_debug,
):
    """Write the file CONFIG and complete the sample names and paths for all
    FASTQ files in PATH.

    PATH is traversed recursively and adds any file with '.fastq' or '.fq' in
    the file name with the file name minus extension as the sample ID.
    """
    if logger_debug:
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)
    
    print_logo()

    # create working dir and db_dir
    os.makedirs(working_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    sample_table = get_samples_from_fastq(path_to_fastq)
    prepare_sample_table_for_atlas(
        sample_table,
        assembler=assembler,
        outfile=os.path.join(working_dir, "samples.tsv"),
    )
    logger.debug(f"\n{sample_table}")

    # Set default binner depending on number of samples
    n_samples = sample_table.shape[0]
    try:
        validate_bingroup_size(sample_table)
    except BinGroupSizeError:
        pass

    create_config(
        db_dir,
        temp_dir,
        os.path.join(working_dir, "config.yaml"),
    )


