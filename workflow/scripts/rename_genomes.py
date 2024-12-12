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

def rename_fasta_file(fasta_in, fasta_out, name_mapping_o2n):
    with open(fasta_in) as ffi, open(fasta_out, "w") as ffo:
        for line in ffi:
           if line[0] == ">":
                for oldvalue, newvalue in name_mapping_o2n:
                    line = line.replace(oldvalue, newvalue)
                ffo.write(f"{line}")
           else:
                ffo.write(line)

def rename_text_file(text_in, text_out, name_mapping_o2n):
    with open(text_in) as tfi, open(text_out, "w") as tfo:
        for line in tfi:
            for oldvalue, newvalue in name_mapping_o2n:
                line = line.replace(oldvalue, newvalue)
            tfo.write(f"{line}")


# Install exception handler
sys.excepthook = handle_exception


# start


from naive_atlas import utils
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
    zip(representatives, utils.gen_names_for_range(len(representatives), prefix=snakemake.params.prefix))
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
mapfile_contigs = snakemake.output.mapfile_contigs
rename_contigs = snakemake.params.rename_contigs


# Bin     Filename        Proteins
paths = pd.read_csv(snakemake.input.paths, sep="\t", index_col=0)


os.makedirs(output_dir)

with open(mapfile_contigs, "w") as out_contigs:
    for rep in representatives:
        new_name = old2new.loc[rep]
        name_mapping_o2n = []
        
        ## Genome
        fasta_in  = paths.loc[rep].Filename
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
                    out_contigs.write(f"{new_header}\t{new_name}\n")
                    name_mapping_o2n.append((line[1:].strip().split()[0], new_header))
                    # write to fasta file
                    ffo.write(f">{new_header}\n")
                else:
                    ffo.write(line)
        
        # Add MAG name to list of renaming pairs - Needs to be last so that we rename scaffolds first, then remaining old MAG names (mostly an issue in GFF)
        name_mapping_o2n.append((rep, new_name))
        
        ## Rename using name mappings
        rename_fasta_file(
            paths.loc[rep].Proteins,
            os.path.join(output_dir, f"{new_name}.faa"),
            name_mapping_o2n,
        )
        rename_fasta_file(
            paths.loc[rep].CDS,
            os.path.join(output_dir, f"{new_name}.fna"),
            name_mapping_o2n,
        )
        if 'rRNA' in paths.columns:
            rename_fasta_file(
                paths.loc[rep].rRNA,
                os.path.join(output_dir, f"{new_name}.rRNA.fna"),
                name_mapping_o2n,
            )
        if 'tRNA' in paths.columns:
            rename_fasta_file(
                paths.loc[rep].tRNA,
                os.path.join(output_dir, f"{new_name}.tRNA.fna"),
                name_mapping_o2n,
            )
        
        rename_text_file(
            paths.loc[rep].GFF,
            os.path.join(output_dir, f"{new_name}.gff"),
            name_mapping_o2n,
        )
        if 'seq_type' in paths.columns:
            rename_text_file(
                paths.loc[rep].seq_type,
                os.path.join(output_dir, f"{new_name}.seq_type.tsv"),
                name_mapping_o2n,
            )


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
