#! /usr/bin/env python

import sys, os, shutil
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


# start
fasta_in       = str(snakemake.input.fa)
fasta_out      = str(snakemake.output.fa)
mapfile_c2g    = str(snakemake.output.mapfile_c2g)
mapfile_g2c    = str(snakemake.output.mapfile_g2c)
rename_contigs = str(snakemake.params.rename_contigs)
new_name       = str(snakemake.params.prefix)

shutil.rmtree(snakemake.params.outdir, ignore_errors=True)
os.makedirs(snakemake.params.outdir)

# write names of contigs in mapping file
with open(fasta_in, "r") as ffi, open(fasta_out, "w") as ffo, open(mapfile_c2g, "w") as mapfile_c2g_fh, open(mapfile_g2c, "w") as mapfile_g2c_fh:
    Nseq = 0
    for line in ffi:
        # if header line
        if line[0] == ">":
            Nseq += 1
            
            if rename_contigs:
                new_header = f"{new_name}-{Nseq:012}"
            else:
                new_header = line[1:].strip().split()[0]
            
            mapfile_c2g_fh.write(f"{new_header}\tUnbinned\n")
            mapfile_g2c_fh.write(f"Unbinned\t{new_header}\n")
            # write to fasta file
            ffo.write(f">{new_header}\n")
        else:
            ffo.write(line)


