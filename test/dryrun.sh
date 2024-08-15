#!/bin/bash
set -euo pipefail


NThreads=72
MaxMem=50

databaseDir="databases"
WD='testrun'
reads_dir='test_reads'
snakemake_args=" --quiet rules --cores $NThreads $@ --dryrun " 
test_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



# Starting clean up
rm -rf $WD $reads_dir $databaseDir "test_reads.tar.gz" ".snakemake" "logs"

echo -e "\n\n\n\n\n\n## Version\n"
naive_atlas --version

echo -e "\n\n\n\n\n\n## Help\n"
naive_atlas run --help

echo -e "\n\n\n\n\n\n## Download read data\n"
wget "https://zenodo.org/record/3992790/files/test_reads.tar.gz"
tar -xzf test_reads.tar.gz

echo -e "\n\n\n\n\n\n## Atlas download\n"
naive_atlas download --db-dir $databaseDir $snakemake_args

echo -e "\n\n\n\n\n\n## Init\n"
naive_atlas init --db-dir $databaseDir --threads=$NThreads -w $WD $reads_dir

echo -e "\n\n\n\n\n\n## Dryrun all\n"
naive_atlas run all -w $WD $snakemake_args



# Ending clean up
#rm -rf $WD $reads_dir $databaseDir "test_reads.tar.gz" ".snakemake" "logs"


