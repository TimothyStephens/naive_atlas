#! /usr/bin/env python

import sys, os
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


from naive_atlas import utils

rename_contigs = snakemake.params.rename_contigs

output_dir = snakemake.output.dir
os.makedirs(output_dir)

for fasta_in in snakemake.input.unbinned:
    new_name = os.path.splitext(os.path.basename(fasta_in))[0]
    
    fasta_out = os.path.join(output_dir, f"{new_name}.fasta")
    
    # write names of contigs in mapping file
    with open(fasta_in) as ffi, open(fasta_out, "w") as ffo:
        Nseq = 0
        for line in ffi:
            # if header line
            if line[0] == ">":
                Nseq += 1
                
                if rename_contigs:
                    new_header = f"{new_name}_{Nseq}"
                else:
                    new_header = line[1:].strip().split()[0]
                
                # write to fasta file
                ffo.write(f">{new_header}\n")
            else:
                ffo.write(line)


