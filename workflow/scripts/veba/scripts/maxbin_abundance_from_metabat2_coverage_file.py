#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import numpy as np
import pandas as pd

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.10.3"


def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <metabat2_coverage_file> -o <output_maxbin_abundance_list_file>".format(__program__)
    epilog = "Copyright 2024 Timothy G. Stephens (ts942@sebs.rutgers.edu)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline

    parser.add_argument("-i","--input",  type=str, help = "path/to/coverage_metabat.tsv [Required]")
    parser.add_argument("-o","--output", type=str, help = "path/to/coverage_maxbin2_list.tsv [Required]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    df = pd.read_csv(opts.input, sep='\t', index_col=0)
    df = df.loc[:,df.columns.map(lambda x: x.endswith('.bam'))]
    
    with open(opts.output, 'w') as paths:
        for i in df.columns.values.tolist():
            print(f"[{__program__}] Processing: {i}")
            p=opts.output+'.PART-'+i
            df[i].to_csv(p, sep='\t', header=None)
            paths.write(p+'\n')
    

if __name__ == "__main__":
    main()
