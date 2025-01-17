#!/usr/bin/env bash

exec 1>"${snakemake_log[0]}" 2>&1 # send all stderr from this script to the log file

mkdir -p "${snakemake_output[dir]}"
for dir in ${snakemake_input[dirs]};
do
    echo "Moving: $dir"
    mv "$dir"/* "${snakemake_output[dir]}"
done
mv -f tmp/genomes/*.genome_quality.tsv "${snakemake_output[dir]}/"
