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


def gen_names_for_range(N, prefix="", start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros = len(str(N))
    format_int = prefix + "{:0" + str(n_leading_zeros) + "d}"
    return [format_int.format(i) for i in range(start, N + start)]


# Install exception handler
sys.excepthook = handle_exception


# start
import pandas as pd

# genome  SpeciesNr       Species Representative
mapping = pd.read_csv(
    snakemake.input.mapping_file,
    sep="\t",
    index_col=0,
).squeeze()


# standardize names of representatives
# MAG001 ....
representatives = mapping.Representative.unique()
old2new_name = dict(
    zip(representatives, gen_names_for_range(len(representatives), prefix=snakemake.params.prefix))
)
mapping["MAG"] = mapping.Representative.map(old2new_name)


# write cluster attribution
mapping[["MAG", "Representative"]].to_csv(
    snakemake.output.mapfile_allbins2mag, sep="\t", header=True
)

# write out old2new ids
old2new = mapping.loc[representatives, "MAG"]
old2new.index.name = "Representative"
old2new.to_csv(snakemake.output.mapfile_old2mag, sep="\t", header=True)


#### Write genomes and contig to genome mapping file
output_dir = snakemake.output.dir
mapfile_c2g = snakemake.output.mapfile_c2g
mapfile_g2c = snakemake.output.mapfile_g2c
rename_contigs = snakemake.params.rename_contigs


# Bin     Genome
paths = pd.read_csv(snakemake.input.paths, sep="\t", index_col=0)


os.makedirs(output_dir)

with open(mapfile_c2g, "w") as mapfile_c2g_fh, open(mapfile_g2c, "w") as mapfile_g2c_fh:
    for rep in representatives:
        new_name = old2new.loc[rep]
        
        ## Genome
        fasta_in  = paths.loc[rep].Genome
        fasta_out = os.path.join(output_dir, f"{new_name}.fa")
        
        # write names of contigs in mapping file
        with open(fasta_in, "r") as ffi, open(fasta_out, "w") as ffo:
            Nseq = 0
            for line in ffi:
                # if header line
                if line[0] == ">":
                    Nseq += 1
                    
                    if rename_contigs:
                        new_header = f"{new_name}-{Nseq:08}"
                    else:
                        new_header = line[1:].strip().split()[0]
                    
                    # write to contig to mapping file
                    mapfile_c2g_fh.write(f"{new_header}\t{new_name}\n")
                    mapfile_g2c_fh.write(f"{new_name}\t{new_header}\n")
                    # write to fasta file
                    ffo.write(f">{new_header}\n")
                else:
                    ffo.write(line)


# rename quality
def rename_quality(quality_in, quality_out, old2new_name):
    Q = pd.read_csv(quality_in, index_col=0, sep="\t")

    Q = Q.loc[old2new_name.keys()].rename(index=old2new_name)

    Q.to_csv(quality_out, sep="\t")


rename_quality(
    quality_in=snakemake.input.genome_info,
    quality_out=snakemake.output.genome_info,
    old2new_name=old2new_name,
)
