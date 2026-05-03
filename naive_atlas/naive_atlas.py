import os
import sys
from .color_logger import logger
from pathlib import Path

import multiprocessing
import subprocess
import click

from snakemake.common.configfile import load_configfile
from .make_config import validate_config
from .init.atlas_init import run_init  # , run_init_sra

from .__init__ import __version__
from .logo import print_logo

##

# Get working dir and scripts dir
cwd = Path.cwd()
script_dir = Path(__file__).resolve().parent

class LogoCommand(click.Command):
    """Click Command subclass that prints the logo before the help message."""
    def format_help(self, ctx, formatter):
        print_logo()
        super().format_help(ctx, formatter)

class LogoGroup(click.Group):
    """Click Group subclass that prints the logo before the help message."""
    def format_help(self, ctx, formatter):
        print_logo()
        super().format_help(ctx, formatter)


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

        # Set Snakemake resources. java_mem is 85% of total to allow for JVM overhead.
        # Total and Java memory are provided in MB for consistency with rule scaling logic.
        total_mb = floor(max_mem * 1024)
        java_mb = floor(0.85 * total_mb)

        return f" --resources 'mem_mb={total_mb}' 'java_mem={java_mb}' "


@click.group(cls=LogoGroup, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """naive_ATLAS - workflows for assembly, annotation, and genomic binning of
    metagenomic and metatranscriptomic data.

    For updates and reporting issues, see: https://github.com/TimothyStephens/naive_atlas
    """


# TODO
# cli.add_command(run_init)
# cli.add_command(run_init_sra)


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


# QC command


@cli.command(
    "run",
    cls=LogoCommand,
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
    )
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas.",
    default=cwd,
)
@click.option(
    "-c",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config-file generated with 'atlas init'",
)
@click.option(
    "--cluster-type",
    type=click.Choice(["slurm", "generic"]),
    default=None,
    help="Type of cluster to use. If set, --profile, --local-cores, and --jobs are sent to Snakemake. Otherwise, --cores is used.",
)
@click.option(
    "--cluster-slurm-params",
    type=str,
    default=f"--slurm-requeue --slurm-efficiency-report --slurm-efficiency-report-path {cwd}{os.sep}efficiency_reports --slurm-efficiency-threshold 1",
    help="Params to use for SLURM cluster. Only used if --cluster-type IS set.",
)
@click.option(
    "--cores",
    type=int,
    default=multiprocessing.cpu_count(),
    help="Number of cores for Snakemake. Defaults to max system cores. Only used if --cluster-type is NOT set.",
)
@click.option(
    "--local-cores",
    type=int,
    default=multiprocessing.cpu_count(),
    help="Number of local cores for Snakemake. Defaults to max system cores. Only used if --cluster-type IS set.",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=250,
    show_default=True,
    help="Use at most this many jobs in parallel. Only used if --cluster-type IS set.",
)
@click.option(
    "--max-mem",
    type=float,
    default=None,
    help=handle_max_mem.__doc__,
)
@click.option(
    "--profile",
    default=None,
    help="Snakemake profile for cluster execution. Default '<cluster-type>_profile'.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.option(
    "--rerun-triggers",
    default="mtime",
    show_default=True,
    help="Define what triggers the rerunning of a job. [{code,input,mtime,params,software-env} ...]",
)
@click.option(
    "--nolock/--lock",
    default=True,
    show_default=True,
    help="Do not lock the working directory.",
)
@click.option(
    "--show-failed-logs/--no-show-failed-logs",
    default=True,
    show_default=True,
    help="Automatically display logs of failed jobs.",
)
@click.option(
    "--scheduler",
    default="greedy",
    show_default=True,
    help="Specifies if jobs are selected by a greedy algorithm or by solving an ilp. [{ilp,greedy}]",
)
@click.option(
    "--printshellcmds/--no-printshellcmds",
    default=True,
    show_default=True,
    help="Print out the shell commands that will be executed.",
)
@click.option(
    "--debug-dag/--no-debug-dag",
    default=True,
    show_default=True,
    help="Print candidate and selected jobs (including their wildcards) while inferring DAG.",
)
@click.option(
    "--verbose/--quiet",
    default=True,
    show_default=True,
    help="Print logging output.",
)
@click.option(
    "--keep-incomplete/--no-keep-incomplete",
    default=True,
    show_default=True,
    help="Do not remove incomplete output files by failed jobs.",
)
@click.option(
    "--keep-going/--no-keep-going",
    default=True,
    show_default=True,
    help="Go on with independent jobs if a job fails.",
)
@click.option(
    "--rerun-incomplete/--no-rerun-incomplete",
    default=True,
    show_default=True,
    help="Re-run jobs that have incomplete output files.",
)
@click.option(
    "--sdm",
    "--software-deployment-method",
    "sdm",
    multiple=True,
    default=["apptainer", "conda"],
    show_default=True,
    help="Software deployment methods (e.g. apptainer, conda).",
)
@click.option(
    "--singularity-args",
    default=f"--no-home --containall --cleanenv --bind {script_dir.parent}:{script_dir.parent} --bind {cwd}:{cwd} ",
    show_default=True,
    help="Arguments passed to singularity/apptainer.",
)
@click.option(
    "--latency-wait",
    type=int,
    default=60,
    show_default=True,
    help="Wait given seconds if an output file of a job is not present after the job finished. This helps if your filesystem suffers from latency.",
)
@click.option(
    "--retries",
    type=int,
    default=1,
    show_default=True,
    help="Number of times to retry failed jobs.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(
    workflow,
    working_dir,
    config_file,
    jobs,
    cores,
    local_cores,
    max_mem,
    profile,
    cluster_type,
    cluster_slurm_params,
    dryrun,
    rerun_triggers,
    nolock,
    show_failed_logs,
    scheduler,
    printshellcmds,
    debug_dag,
    verbose,
    keep_incomplete,
    keep_going,
    rerun_incomplete,
    sdm,
    singularity_args,
    latency_wait,
    retries,
    snakemake_args,
):
    """Runs the naive ATLAS pipline
    
    By default all steps are executed but a sub-workflow can be specified.
    Needs a config-file and expects to find a sample table in the working-directory. Both can be generated with 'atlas init'
    
    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    
    \b
    # OPTIONS:
    | qc -> assembly -> binning -> genomes -> quantify_genomes ---------------------------------------|-> strains
    |                                    +-> genome_annotation -------------------------------------->|
    |                                                        +-> gene_prediction -> gene_annotation ->|
    +-----------------------------------------------all-----------------------------------------------+
    # Independent of other steps:
    screen
    download (download reference databases upfront instead of as each rule needs them (need ~920GB for all databases, ~1.5TB during download))

    """
    
    print_logo()
    logger.info("STARTING WORKFLOW!")

    cluster_params = ""
    if cluster_type:
        if profile is None:
            profile = cluster_type + '_profile'
        core_str = " --local-cores {} --jobs {} ".format(
            local_cores, jobs
        )
        if cluster_type == "slurm":
            cluster_params = cluster_slurm_params
    else:
        core_str = " --cores {} ".format(cores)

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

    validate_config(config_file, workflow)

    # Helper to generate presence-only flags
    def get_flag(val, name):
        return f" --{name} " if val else ""

    cmd = (
        "snakemake "

        # Snakemake setup
        " --snakefile '{snakefile}' "
        " --configfile '{config_file}' "
        " --directory '{working_dir}' "
        " {profile} "
        " {target_rule} "
        
        # Snakemake run behavior params
        " --rerun-triggers {rerun_triggers} "
        " {nolock} "
        " --scheduler {scheduler} "
        " {keep_going} "
        " {rerun_incomplete} "
        " {keep_incomplete} "
        " --latency-wait {latency} "
        " --retries {retries} "
        
        # Logging
        " {printshellcmds} "
        " {show_failed_logs} "
        " {verbose} "
        " {debug_dag} "

        # Env setup
        " --software-deployment-method {sdm} "
        " --singularity-args '{sing_args}' "

        # Resource limits
        " {core_str} "
        " {max_mem_string} "

        # Extra params & dryrun
        " {cluster_params} "
        " {args} "
        " {dryrun} "
    ).format(
        # Snakemake setup
        snakefile=get_snakefile(),
        config_file=config_file,
        working_dir=working_dir,
        profile="" if (profile is None) else "--profile {}".format(profile),
        target_rule=workflow if workflow != "None" else "",

        # Snakemake run behavior params
        rerun_triggers=rerun_triggers,
        nolock=get_flag(nolock, "nolock"),
        scheduler=scheduler,
        keep_going=get_flag(keep_going, "keep-going"),
        rerun_incomplete=get_flag(rerun_incomplete, "rerun-incomplete"),
        keep_incomplete=get_flag(keep_incomplete, "keep-incomplete"),
        latency=latency_wait,
        retries=retries,

        # Logging
        printshellcmds=get_flag(printshellcmds, "printshellcmds"),
        show_failed_logs=get_flag(show_failed_logs, "show-failed-logs"),
        verbose=get_flag(verbose, "verbose"),
        debug_dag=get_flag(debug_dag, "debug-dag"),
        
        # Env setup
        sdm=" ".join(sdm),
        sing_args=singularity_args,

        # Resource limits
        core_str=f"{core_str}",
        max_mem_string=handle_max_mem(max_mem, profile),

        # Extra params & dryrun
        cluster_params=f"{cluster_params}",
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
    
    logger.info("FINISHED WORKFLOW!")

if __name__ == "__main__":
    cli()
