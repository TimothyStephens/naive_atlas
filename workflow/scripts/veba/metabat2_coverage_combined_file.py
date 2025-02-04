#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd
from functools import reduce

pd.options.display.max_colwidth = 100
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.01.23"


def main(args=None):
    # Program Info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <coverage_file_1> <coverage_file_2> ... -o <output_metabat2_coverage_file>".format(__program__)
    epilog = "Copyright 2025 Timothy G. Stephens (ts942@sebs.rutgers.edu)"
    
    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--input",
            required=True, type=str, nargs='+',
            help = "path/to/coverage_file_1.tsv path/to/coverage_file_2.tsv ... [Required]"
    )
    parser.add_argument("-o","--output",
            required=True, type=str,
            help = "path/to/coverage_metabat2.tsv [Required]"
    )
    
    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Load each depth file and merge using columns which are identical across files: 'contigName', 'contigLen', 'totalAvgDepth'
    data_frames = [ pd.read_table(file_name, sep='\t') for file_name in opts.input ]
    df = data_frames[0]
    for d in data_frames[1:]:
        # totalAvgDepth is different in each file and needs to be recaluclated at the end
        d.drop('totalAvgDepth', axis=1, inplace=True)
        df = df.merge(d, on=['contigName', 'contigLen'], how='outer')
    
    # Calculate the average of NON-"var" columns
    cols_to_average = [col for col in df.columns if not col.endswith('-var') and not col in ['contigName', 'contigLen', 'totalAvgDepth'] ]
    df['totalAvgDepth'] = df[cols_to_average].mean(axis=1)
    
    df.to_csv(opts.output, sep='\t', index=False, na_rep=0.0)


if __name__ == "__main__":
    main()
