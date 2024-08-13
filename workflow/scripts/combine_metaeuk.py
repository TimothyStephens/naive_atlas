import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


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

import pandas as pd


def main(genomes, contig_fasta_files, contig_results_files, mag_results_files, contig_output_table, mag_output_table):
    ## Contig
    df_list = []
    for i, genome in enumerate(genomes):
        headers = [x.strip().lstrip('>').split(' ')[0] for x in open(contig_fasta_files[i], 'r').readlines() if x.startswith('>')]
        
        all_names = pd.DataFrame({'contig_names':headers}, index=headers)
        
        data = pd.read_table(contig_results_files[i], index_col=0, dtype=str)
        data = pd.concat([all_names, data], axis=1)
        data["Genome"] = genome
        
        df_list.append(data)

    df = pd.concat(df_list, axis=0)
    
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    
    df = df.fillna('unclassifiable')
    
    df.to_csv(contig_output_table, sep="\t", index=False)

    ## MAG
    df_list = []
    for i, genome in enumerate(genomes):
        data = pd.read_table(mag_results_files[i], index_col=False, dtype=str)

        df_list.append(data)

    df = pd.concat(df_list, axis=0)

    df = df.fillna('unclassifiable')

    df.to_csv(mag_output_table, sep="\t", index=False)



if __name__ == "__main__":
    main(
        genomes=snakemake.params.genomes,
        contig_fasta_files=snakemake.input.contig_fasta_files,
        contig_results_files=snakemake.input.contig_results_files,
        mag_results_files=snakemake.input.mag_results_files,
        contig_output_table=snakemake.output.contig_output_table,
        mag_output_table=snakemake.output.mag_output_table,
    )
