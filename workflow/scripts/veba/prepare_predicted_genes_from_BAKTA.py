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
    usage = "{} -i <BAKTA/PROKKA_pep_1> <BAKTA/PROKKA_pep_2> -o <output_directory> -s <stats_file>".format(__program__)
    epilog = "Copyright 2025 Timothy G. Stephens"
    
    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Pipeline
    parser.add_argument("-i","--input_pep_files",  type=str, required=True, help="path/to/bakta_output_files/*.faa", nargs='*')
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
        in_faa = open_file(prefix + ".faa", 'r')
        in_fna = open_file(prefix + ".ffn", 'r')
        in_gff = open_file(prefix + ".gff3", 'r')
        in_tsv = open_file(prefix + ".tsv", 'r')
        in_txt = open_file(prefix + ".txt", 'r')
        
        # Output files
        out_faa       = open_file(os.path.join(opts.output_directory, name + ".faa"),  'w')
        out_fna       = open_file(os.path.join(opts.output_directory, name + ".fna"),  'w')
        out_gff       = open_file(os.path.join(opts.output_directory, name + ".gff3"), 'w')
        out_tsv       = open_file(os.path.join(opts.output_directory, name + ".tsv"),  'w')
        out_rRNA_fna  = open_file(os.path.join(opts.output_directory, name + ".rRNA.fna"),   'w')
        out_rRNA_gff  = open_file(os.path.join(opts.output_directory, name + ".rRNA.gff3"),  'w')
        out_tRNA_fna  = open_file(os.path.join(opts.output_directory, name + ".tRNA.fna"),   'w')
        out_tRNA_gff  = open_file(os.path.join(opts.output_directory, name + ".tRNA.gff3"),  'w')
        out_ncRNA_fna = open_file(os.path.join(opts.output_directory, name + ".ncRNA.fna"),  'w')
        out_ncRNA_gff = open_file(os.path.join(opts.output_directory, name + ".ncRNA.gff3"), 'w')
        out_other_fna = open_file(os.path.join(opts.output_directory, name + ".other.fna"),  'w')
        out_other_gff = open_file(os.path.join(opts.output_directory, name + ".other.gff3"), 'w')
        
        # Parse GFF to extract header
        gff_header = []
        for line in in_gff:
            line = line.rstrip('\n')
            if line.startswith('##sequence-region'):
                break
            gff_header.append(line)
        in_gff.seek(0)
        
        out_gff.write('\n'.join(gff_header) + '\n')
        out_rRNA_gff.write('\n'.join(gff_header) + '\n')
        out_tRNA_gff.write('\n'.join(gff_header) + '\n')
        out_ncRNA_gff.write('\n'.join(gff_header) + '\n')
        out_other_gff.write('\n'.join(gff_header) + '\n')
        
        # Parse GFF to extract different feature types
        rRNA_names  = []
        tRNA_names  = []
        ncRNA_names = []
        other_names = []
        for line in in_gff:
            line = line.rstrip('\n')
            line_split = line.split('\t')
            
            if len(line_split) != 9 or line.startswith('#') or line_split[2] == "region":
                continue
            
            features = { f.split('=')[0]:f.split('=')[1] for f in line_split[8].split(';') }
            
            if line_split[2] == "CDS":
                out_gff.write(f'{line}\n')
            
            elif line_split[2] == "rRNA":
                out_rRNA_gff.write(f'{line}\n')
                rRNA_names.append(features['locus_tag'])
            
            elif line_split[2] == "tRNA":
                out_tRNA_gff.write(f'{line}\n')
                tRNA_names.append(features['locus_tag'])
            
            elif line_split[2] == "ncRNA":
                out_ncRNA_gff.write(f'{line}\n')
                ncRNA_names.append(features['locus_tag'])
            
            else: # tmRNA, CRISPR, etc.
                out_other_gff.write(f'{line}\n')
                other_names.append(features['ID'])
        
        # Parse CDS to extract different feature types
        for header, seq in SimpleFastaParser(in_fna):
            header = header.split(' ')[0]
            
            if header in rRNA_names:
                out_rRNA_fna.write(f'>{header}\n{seq}\n')
            
            elif header in tRNA_names:
                out_tRNA_fna.write(f'>{header}\n{seq}\n')
            
            elif header in ncRNA_names:
                out_ncRNA_fna.write(f'>{header}\n{seq}\n')
            
            elif header in other_names:
                out_other_fna.write(f'>{header}\n{seq}\n')
            
            else:
                out_fna.write(f'>{header}\n{seq}\n')
        
        # Parse PEP to extract different feature types
        for header, seq in SimpleFastaParser(in_faa):
            header = header.split(' ')[0]
            out_faa.write(f'>{header}\n{seq}\n')
        
        # Parse gene description TSV
        for line in in_tsv:
            line = line.strip('\n')
            if line.startswith('#Sequence'):
                line = line.lstrip('#')
            out_tsv.write(f'{line}\n')
        
        # Parse gene stats TXT
        stats = {}
        stats["MAG"] = name
        for line in in_txt:
            line = line.strip('\n')
            if not line or line.startswith('#'):
                continue
            
            line_split = line.split(': ')
            if len(line_split) == 2:
                stats[line_split[0]] = line_split[1]
        del stats['Software']
        del stats['Database']
        del stats['DOI']
        del stats['URL']
        all_stats.append(stats)
        
        # Close everything
        in_faa.close()
        in_fna.close()
        in_gff.close()
        in_tsv.close()
        in_txt.close()
        out_faa.close()
        out_fna.close()
        out_gff.close()
        out_tsv.close()
        out_rRNA_fna.close()
        out_rRNA_gff.close()
        out_tRNA_fna.close()
        out_tRNA_gff.close()
        out_ncRNA_fna.close()
        out_ncRNA_gff.close()
        out_other_fna.close()
        out_other_gff.close()
    
    df = pd.DataFrame(all_stats)
    df.set_index('MAG', inplace=True)
    df.to_csv(opts.stats_file, sep='\t', quoting=False)

if __name__ == "__main__":
    main()                
