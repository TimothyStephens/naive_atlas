import os, sys, shutil
from .color_logger import logger

import multiprocessing
import subprocess
import click

from snakemake.common.configfile import load_configfile
from .__init__ import __version__



###### Functions ######
def handle_max_mem(max_mem, profile):
    "Specify maximum virtual memory to use by atlas."
    "For numbers >1 its the memory in GB. "
    "For numbers <1 it's the fraction of available memory."

    if profile is not None:
        if max_mem is not None:
            logger.info(
                "Memory requirements are handled by the profile, I ignore max-mem argument."
            )
        # memory is handled via the profile, user should know what he is doing
        return ""
    else:
        import psutil
        from math import floor

        # calulate max  system meory in GB (float!) - Include SWAP for systems that have it avaliale
        max_system_memory = (psutil.virtual_memory().total / (1024**3)) + (psutil.swap_memory().total / (1024**3))

        if max_mem is None:
            max_mem = 0.95
        if max_mem > 1:
            if max_mem > max_system_memory:
                logger.critical(
                    f"You specified {max_mem} GB as maximum memory, but your system only has {floor(max_system_memory)} GB"
                )
                sys.exit(1)

        else:
            max_mem = max_mem * max_system_memory

        # specify max_mem_string including java mem and max mem

        return f" --resources mem={floor(max_mem)} mem_mb={floor(max_mem*1024)} java_mem={floor(0.85* max_mem)} "


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf



###### Naive ATLAS ######

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """naive_ATLAS - workflows for assembly, annotation, and genomic binning of
    metagenomic and metatranscriptomic data.

    For updates and reporting issues, see: https://github.com/TimothyStephens/naive_atlas
    """


###### init command ######


@cli.command(
    "init",
    short_help="prepare configuration file and sample table for atlas run",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas",
    default=".",
)
def run_init(
    working_dir,
):
    """Write the CONFIG and example SAMPLES files to working dir.
    """
    config_file  = os.path.join(working_dir, "config.yaml")
    samples_file = os.path.join(working_dir, "samples.tsv")

    # Create working dir
    os.makedirs(working_dir, exist_ok=True)

    # Create template config.yaml file
    shutil.copy2(
        "../config/template_samples.tsv", # Relative from naive_atlas/naive_atlas/naive_atlas.py
        samples_file
    )
    print()
    print(f"## Created template config file: {config_file}")
    print(f"## Please add host datasets for mapping (if any exist) and uncomment the annotation approaches that you want to run (it is also fine to keep just the default config and change nothing).")
    print()

    # Create template samples.tsv file
    with open("../config/template_config.yaml", 'r') as infile: # Relative from naive_atlas/naive_atlas/naive_atlas.py
        file_content = infile.read()
    new_content = file_content.replace('/project', working_dir)
    with open(config_file, 'w') as outfile:
        outfile.write(new_content)
    print()
    print(f"## Created template samples file: {samples_file}")
    print(f"## Please use this template to guide the construction of your final samples file.")
    print()



###### run command ######
@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run atlas main workflow",
)
@click.argument(
    "workflow",
    type=click.Choice(
        [
            "download",
            "qc",
            "assembly",
            "binning",
            "genomes",
            "quantify_genomes",
            "genome_annotation",
            "gene_prediction",
            "gene_annotation",
            "strains",
            "screen",
            "None",
            "all",
            "test",
        ]
    ),
    #    show_default=True,
    #    help="Execute only subworkflow.",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas.",
    default=".",
)
@click.option(
    "-c",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config-file generated with 'atlas init'",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for more details).",
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "--max-mem",
    type=float,
    default=None,
    help=handle_max_mem.__doc__,
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(
    workflow, working_dir, config_file, jobs, max_mem, profile, dryrun, snakemake_args
):
    """Runs the naive_ATLAS pipline
    
    By default all steps are executed but a sub-workflow can be specified.
    Needs a config-file and expects to find a sample table in the working-directory. Both can be generated with 'atlas init'
    
    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    
    For more details, see: https://metagenome-atlas.readthedocs.io
    
    \b
    # OPTIONS:
    | qc -> assembly -> binning -> genomes -> quantify_genomes ---------------------------------------|-> strains
    |                                    +-> genome_annotation -------------------------------------->|
    |                                                        +-> gene_prediction -> gene_annotation ->|
    +-----------------------------------------------all-----------------------------------------------+
    # Independent of other steps:
    screen
    download (download reference files (need ~920GB for all databases, ~1.5TB during download))

    """

    logger.info(f"Atlas version: {__version__}")

    if config_file is None:
        config_file = os.path.join(working_dir, "config.yaml")

    if not os.path.exists(config_file):
        logger.critical(
            f"config-file not found: {config_file}\n" "generate one with 'atlas init'"
        )
        exit(1)

    sample_file = os.path.join(working_dir, "samples.tsv")

    if not os.path.exists(sample_file):
        logger.critical(
            f"sample.tsv not found in the working directory. "
            "Generate one with 'atlas init'"
        )
        exit(1)

    conf = load_configfile(config_file)

    db_dir = conf["database_dir"]

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        " --configfile '{config_file}' "
        " {profile} "
        "{jobs} {max_mem_string} "
        " --keep-incomplete --rerun-incomplete --keep-going "
        " --rerun-triggers mtime "
        " --show-failed-logs "
        " --use-conda --use-apptainer "
        " --scheduler greedy "
        " {target_rule} "
        " {args} "
        " {dryrun} "
    ).format(
        snakefile=get_snakefile(),
        working_dir=working_dir,
        config_file=config_file,
        profile="" if (profile is None) else "--profile {}".format(profile),
        jobs="--jobs {}".format(jobs) if jobs is not None else "",
        max_mem_string=handle_max_mem(max_mem, profile),
        target_rule=workflow if workflow != "None" else "",
        args=" ".join(snakemake_args),
        dryrun="--dryrun" if dryrun else "",
    )
    logger.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)


if __name__ == "__main__":
    cli()
