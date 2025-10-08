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


timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%X")


def get_read_stats(fraction, params_in):
    "get read stats by running reformat.sh"

    from snakemake.shell import shell

    subfolder = os.path.join(snakemake.params.folder, fraction)
    tmp_file = os.path.join(subfolder, "read_stats.tmp")
    ## `qhist` is commented out becuase it can cause the following error with single-end reads, for some unknown reason.
    # Exception in thread "main" java.lang.AssertionError: NaN, 0.0, 1.0
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # tracker.ReadStats@7cca494b
    #	at tracker.ReadStats.calcEntropySuperSlow(ReadStats.java:1316)
    #	at tracker.ReadStats.writeQualityToFile(ReadStats.java:1046)
    #	at tracker.ReadStats.writeAll(ReadStats.java:852)
    #	at jgi.ReformatReads.process(ReformatReads.java:1210)
    #	at jgi.ReformatReads.main(ReformatReads.java:54)
    shell(
        f" mkdir -p {subfolder} 2>> {snakemake.log[0]} "
        f" ; "
        f" reformat.sh {params_in} "
        f" bhist={subfolder}/base_hist.txt "
        #f" qhist={subfolder}/quality_by_pos.txt " 
        f" lhist={subfolder}/readlength.txt "
        f" gchist={subfolder}/gc_hist.txt "
        f" gcbins=auto "
        f" bqhist={subfolder}/boxplot_quality.txt "
        f" threads={snakemake.threads} "
        f" overwrite=true "
        f" -Xmx{snakemake.resources.java_mem}G "
        f" 2>&1 | tee -a {snakemake.log[0]} {tmp_file} >/dev/null "
    )
    content = open(tmp_file).read()
    pos = content.find("Input:")
    if pos == -1:
        raise Exception("Didn't find read number in file:\n\n" + content)
    else:
        content[pos:].split()[1:4]
        # Input:    123 reads   1234 bases
        n_reads, _, n_bases = content[pos:].split()[1:4]

        os.remove(tmp_file)
    return int(n_reads), int(n_bases)


# Generate stats for each read type if present.
if hasattr(snakemake.input, 'R1'):
    n_reads_pe, n_bases_pe = get_read_stats(
        "pe", f"in1={snakemake.input.R1} in2={snakemake.input.R2}"
    )
    n_reads_pe = n_reads_pe / 2
else:
    n_reads_pe, n_bases_pe = 0, 0


if hasattr(snakemake.input, 'SE'):
    n_reads_se, n_bases_se = get_read_stats(
        "se", f"in={snakemake.input.SE}"
    )
else:
    n_reads_se, n_bases_se = 0, 0


if hasattr(snakemake.input, 'LR'):
    n_reads_lr, n_bases_lr = get_read_stats(
        "lr", f"in={snakemake.input.LR}"
    )
else:
    n_reads_lr, n_bases_lr = 0, 0


headers = [
    "Sample",
    "Step",
    "Total_Reads",
    "Total_Bases",
    "Reads_pe",
    "Bases_pe",
    "Reads_se",
    "Bases_se",
    "Reads_lr",
    "Bases_lr",
    "Timestamp",
]

values = [
    n_reads_pe + n_reads_se + n_reads_lr,
    n_bases_pe + n_bases_se + n_bases_lr,
    n_reads_pe,
    n_bases_pe,
    n_reads_se,
    n_bases_se,
    n_reads_lr,
    n_bases_lr,
]


with open(snakemake.output.read_counts, "w") as f:
    f.write("\t".join(headers) + "\n")
    f.write(
        "\t".join(
            [snakemake.wildcards.sample, snakemake.wildcards.step]
            + [str(v) for v in values]
            + [timestamp]
        )
        + "\n"
    )

shutil.make_archive(snakemake.params.folder, "zip", snakemake.params.folder)
shutil.rmtree(snakemake.params.folder)
