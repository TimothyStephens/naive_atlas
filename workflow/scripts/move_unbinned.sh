#!/usr/bin/env bash

exec 1>"${snakemake_log[0]}" 2>&1 # send all stderr from this script to the log file

mkdir -p "${snakemake_output[dir]}"
for fa in ${snakemake_input[fa]};
do
    echo "Moving: $fa"
    mv "$fa" "${snakemake_output[dir]}"
done
