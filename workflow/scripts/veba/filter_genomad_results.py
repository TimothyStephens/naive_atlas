#!/usr/bin/env python
import sys, os, glob, argparse
from shutil import copyfile, rmtree
# from collections import OrderedDict
import pandas as pd
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.2.14"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <checkv_output.tsv> -f <scaffolds.fasta> -o <output_directory>  -m 1000".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("--scaffolds_to_bins", type=str, required=True, help = "path/to/binning/scaffolds_to_bins.tsv")
    parser.add_argument("-i","--genomad_results", type=str, required=True, help = "path/to/checkv/{virus,plasmid}_summary.tsv")
    parser.add_argument("--genomad_virus_taxonomy", type=str, help = "geNomad virus_taxonomy.tsv")
    parser.add_argument("-f", "--fasta", required=True, type=str, help = "path/to/scaffolds.fasta")
    parser.add_argument("--output_directory", type=str, default="bins_filtered", help = "path/to/bins_filtered [Default: bins_filtered]")
    parser.add_argument("-p", "--prefix", type=str, default="Virus.", help = "Prefix for viral name (e.g. Sample_1__Virus. would make Sample_1__Virus.1, Sample_2__Virus.2, etc.) [Default: Virus.]")
    parser.add_argument("-u", "--unbinned", action="store_true", help="Write unbinned fasta sequences to file")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Output directory
    rmtree(opts.output_directory, ignore_errors=True)
    os.makedirs(opts.output_directory, exist_ok=True)
    os.makedirs(os.path.join(opts.output_directory,"genomes"), exist_ok=True)

    # Load genomad
    df_genomad = pd.read_csv(opts.genomad_results, sep="\t", index_col=0)
    # Add taxonomy if avaliable
    if not opts.genomad_virus_taxonomy is None:
        df_genomad.drop('taxonomy', axis=1, inplace=True)
        df_genomad_virus_taxonomy = pd.read_csv(opts.genomad_virus_taxonomy, sep="\t", index_col=0)
        df_genomad = pd.merge(df_genomad, df_genomad_virus_taxonomy, on='seq_name', how='left')
    
    # Filter functions
    if not opts.genomad_virus_taxonomy is None:
        def filter_genomad(series):
            length, n_genes, virus_score, fdr, n_hallmarks, marker_enrichment, n_genes_with_taxonomy, agreement, lineage = series[["length", "n_genes", "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "n_genes_with_taxonomy", "agreement", "lineage"]]
            conditions = [
                length > 0,
                n_genes > 0,
                virus_score > 0,
                fdr < 0.05,
                n_hallmarks > 0,
                marker_enrichment > 0,
                n_genes_with_taxonomy > 0,
                agreement > 0,
                not lineage is None,
            ]

            return all(conditions)
    else:
        def filter_genomad(series):
            length, n_genes, plasmid_score, fdr, n_hallmarks, marker_enrichment = series[["length", "n_genes", "plasmid_score", "fdr", "n_hallmarks", "marker_enrichment"]]
            conditions = [
                length > 0,
                n_genes > 0,
                plasmid_score > 0,
                fdr < 0.05,
                n_hallmarks > 0,
                marker_enrichment > 0,
            ]

            return all(conditions)

    # Filtered Results
    mask = df_genomad.apply(filter_genomad, axis=1)
    df_genomad = df_genomad.loc[mask]
    
    # Sort
    df_genomad = df_genomad.sort_values('length', ascending=False)
    # Add MAG ID
    df_genomad['file'] = df_genomad.reset_index().index+1
    df_genomad['file'] = opts.prefix+df_genomad['file'].astype(str)
    # Reorder
    cols = list(df_genomad)
    cols.insert(0, cols.pop(cols.index('file')))
    df_genomad = df_genomad.loc[:, cols]
    
    # Contig lists
    f_binned_list = open(os.path.join(opts.output_directory, "binned.list"), "w")
    f_unbinned_list = open(os.path.join(opts.output_directory, "unbinned.list"), "w")
    if opts.unbinned:
        f_unbinned_fasta = open(os.path.join(opts.output_directory, "unbinned.fasta"), "w")
    else:
        f_unbinned_fasta = open(os.devnull, "w")
    if df_genomad.empty:
        print("No MAGs remain after filtering.")
        with open(opts.fasta, "r") as f_fasta: # Use stdin?
            for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting viral and unbinned contigs", unit=" contig"):
                id_scaffold = header.split(" ")[0]
                print(id_scaffold, file=f_unbinned_list)
                print(">{}\n{}".format(header, seq), file=f_unbinned_fasta)
        df_genomad = pd.DataFrame(columns=["id_contig"] + df_genomad.columns.tolist())

    else:
        print("MAGs remain after filtering.")
        # Quality assessment on MAGs
        scaffold_to_bin = dict()
        with open(opts.scaffolds_to_bins, "r") as f_s2b:
            for l in f_s2b:
                if not l:
                    continue
                s, b = l.strip().split('\t')
                if b in df_genomad.index:
                    scaffold_to_bin[s] = b
        
        print(scaffold_to_bin)
        with open(opts.fasta, "r") as f_fasta: # Use stdin?
            for header, seq in tqdm(SimpleFastaParser(f_fasta), "Extracting viral and unbinned contigs", unit=" contig"):
                id_scaffold = header.split(" ")[0]

                if id_scaffold in scaffold_to_bin:
                    id_name = "{}".format(df_genomad.loc[scaffold_to_bin[id_scaffold], "file"])
                    with open(os.path.join(opts.output_directory, "genomes", "{}.fa".format(id_name)), "a") as f_out:
                        scaffold_to_bin[id_scaffold] = id_name
                        print(id_scaffold, file=f_binned_list)
                        print(">{}\n{}".format(id_scaffold, seq), file=f_out)
                else:
                    print(id_scaffold, file=f_unbinned_list)
                    print(">{}\n{}".format(header, seq), file=f_unbinned_fasta)
        scaffold_to_bin = pd.Series(scaffold_to_bin)
        scaffold_to_bin.to_frame().to_csv(os.path.join(opts.output_directory, "scaffolds_to_bins.tsv"), sep="\t", header=None)

        with open(os.path.join(opts.output_directory, "bins.list"), "w") as f_bins:
            for id_name in scaffold_to_bin.unique():
                print(id_name, file=f_bins)

    # Output table
    df_output = df_genomad
    
    df_output.index.name = "file"
    df_output.to_csv(os.path.join(opts.output_directory, "genomad_results.filtered.tsv" ), sep="\t", index=False, na_rep='NA')
    
    f_binned_list.close()
    f_unbinned_list.close()
    f_unbinned_fasta.close()

if __name__ == "__main__":
    main()
    
                
