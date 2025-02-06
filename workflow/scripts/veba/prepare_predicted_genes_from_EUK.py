#!/usr/bin/env python
import sys, os, glob, argparse, gzip
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.2.5"

def open_file(file_name, mode):
    print(f'Opening {file_name} in mode {mode}')
    
    ## Check file exists if mode='r'
    if not os.path.exists(file_name) and mode == 'r':
        raise Exception("The file %s does not exist!" % file_name)
    
    ## Open with gzip if it has the *.gz extension, else open normally (including stdin)
    try:
        if file_name.endswith(".gz"):
            #print "Opening gzip compressed file (mode: %s): %s" % (mode, file_name) ## DEBUG
            file_obj = gzip.open(file_name, mode+'b')
        else:
            #print "Opening normal file (mode: %s): %s" % (mode, file_name) ## DEBUG
            file_obj = open(file_name, mode)
    except IOError as e:
        raise Exception('%s' % e)
    
    return(file_obj)


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <VEBA_pep_1> <VEBA_pep_2> -o <output_directory> -s <stats_file>".format(__program__)
    epilog = "Copyright 2025 Timothy G. Stephens"
    
    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Pipeline
    parser.add_argument("-i","--input_pep_files",  type=str, required=True, help="path/to/VEBA_output_files/*.faa", nargs='*')
    parser.add_argument("-o","--output_directory", type=str, required=True, help="path/to/output_directory/")
    parser.add_argument("-s","--stats_file",       type=str, required=True, help="path/to/stats_file.tsv")
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Output dir (if it doesnt exist)
    os.makedirs(opts.output_directory, exist_ok=True)
    
    file_names = opts.input_pep_files
    file_names = [f.rstrip('.faa') for f in file_names]
    print(f'Found {len(file_names)} gene sets to process.')
    if len(file_names) == 0:
        print("Stopping early. No files to process.")
        return()
    
    all_stats = []
    for prefix in file_names:
        name = os.path.basename(prefix)
        print(f'Processing: {prefix}')
        
        # Input files
        in_faa  = open_file(prefix + ".faa", 'r')
        in_fna  = open_file(prefix + ".fna", 'r')
        in_gff  = open_file(prefix + ".gff", 'r')
        in_rRNA = open_file(prefix + ".rRNA.fna", 'r')
        in_tRNA = open_file(prefix + ".tRNA.fna", 'r')
        in_seq  = open_file(prefix + ".seq_type.tsv", 'r')
        
        # Output files
        out_faa       = open_file(os.path.join(opts.output_directory, name + ".faa"),  'w')
        out_fna       = open_file(os.path.join(opts.output_directory, name + ".fna"),  'w')
        out_gff       = open_file(os.path.join(opts.output_directory, name + ".gff3"), 'w')
        out_rRNA_fna  = open_file(os.path.join(opts.output_directory, name + ".rRNA.fna"),     'w')
        out_rRNA_gff  = open_file(os.path.join(opts.output_directory, name + ".rRNA.gff3"),    'w')
        out_tRNA_fna  = open_file(os.path.join(opts.output_directory, name + ".tRNA.fna"),     'w')
        out_tRNA_gff  = open_file(os.path.join(opts.output_directory, name + ".tRNA.gff3"),    'w')
        out_seq       = open_file(os.path.join(opts.output_directory, name + ".seq_type.tsv"), 'w')
        
        # Parse GFF to extract header
        gff_header = []
        for line in in_gff:
            line = line.rstrip('\n')
            if not line.startswith('#'):
                break
            gff_header.append(line)
        in_gff.seek(0)
        
        out_gff.write('\n'.join(gff_header) + '\n')
        out_rRNA_gff.write('\n'.join(gff_header) + '\n')
        out_tRNA_gff.write('\n'.join(gff_header) + '\n')
        
        gene_count = 0
        exon_count = 0
        rRNA_count = 0
        tRNA_count = 0
        # Parse GFF to extract different feature types
        for line in in_gff:
            line = line.rstrip('\n')
            line_split = line.split('\t')
            
            if len(line_split) != 9 or line.startswith('#') or line_split[2] == "region":
                continue
            
            if line_split[1] == "MetaEuk":
                if line_split[2] == "gene":
                    gene_count+=1
                    out_gff.write(f'{line}\n')
                if line_split[2] == "exon":
                    exon_count+=1
                    # Reverse coords if negative strand since MetaEuk output this way
                    if line_split[6] == '-':
                        t = line_split[3]
                        line_split[3] = line_split[4]
                        line_split[4] = t
                    t = "\t".join(line_split)
                    line_split[2] = "CDS"
                    t = "\t".join(line_split)+ '\n' + t
                    out_gff.write(f'{t}\n')
            
            elif line_split[1].startswith("pyrodigal"):
                    gene_count+=1
                    exon_count+=1
                    out_gff.write(f'{line}\n')
            
            elif line_split[1].startswith("barrnap"):
                if line_split[2] == "rRNA":
                    rRNA_count+=1
                out_rRNA_gff.write(f'{line}\n')
            
            elif line_split[1] == "tRNAscan-SE":
                if line_split[2] == "tRNA":
                    tRNA_count+=1
                out_tRNA_gff.write(f'{line}\n')
            
            else: # Should not hit, warn user if we do
                print(f'[WARNING] Found unexpected gff feature: {line}')
                out_gff.write(f'{line}\n')
        
        # Parse CDS to extract different feature types
        for header, seq in SimpleFastaParser(in_fna):
            header = header.split(' ')[0]
            out_fna.write(f'>{header}\n{seq}\n')
        
        # Parse PEP to extract different feature types
        for header, seq in SimpleFastaParser(in_faa):
            header = header.split(' ')[0]
            out_faa.write(f'>{header}\n{seq}\n')
        
        # Parse seq description TSV
        out_seq.write(f'scaffold\tseq_type\n')
        for line in in_seq:
            line = line.strip('\n')
            if line.startswith('#Sequence'):
                line = line.lstrip('#')
            out_seq.write(f'{line}\n')
        
        # Parse gene stats TXT
        stats = {}
        stats["MAG"] = name
        stats["genes"] = gene_count
        stats["exons"] = exon_count
        stats["rRNAs"] = rRNA_count
        stats["tRNAs"] = tRNA_count
        all_stats.append(stats)
        
        # Close everything
        in_faa.close()
        in_fna.close()
        in_gff.close()
        in_seq.close()
        out_faa.close()
        out_fna.close()
        out_gff.close()
        out_rRNA_fna.close()
        out_rRNA_gff.close()
        out_tRNA_fna.close()
        out_tRNA_gff.close()
        out_seq.close()
    
    df = pd.DataFrame(all_stats)
    df.set_index('MAG', inplace=True)
    df.to_csv(opts.stats_file, sep='\t', quoting=False)

if __name__ == "__main__":
    main()                
