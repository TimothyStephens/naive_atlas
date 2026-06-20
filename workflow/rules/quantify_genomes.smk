

#####################################
####                             ####
#### Genome mapping and coverage ####
####                             ####
#####################################

def get_all_genomes(wildcards):
    # check if genomes are present
    genomes = glob_wildcards("genomes/genomes/{genome}.fa").genome
    
    if len(genomes) == 0 and os.path.isdir("genomes/genomes"):
        logger.error(
            f"No genomes found with fasta extension in genomes/genomes "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        #exit(1)

    return genomes


def get_all_unbinned(wildcards):
    # check if genomes are present
    genomes = glob_wildcards("genomes/unbinned/{genome}.fa").genome

    if len(genomes) == 0 and os.path.isdir("genomes/unbinned"):
        logger.error(
            f"No genomes found with fasta extension in genomes/unbinned "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        #exit(1)

    return genomes


### Quantification
localrules:
    concat_genomes,

rule concat_genomes:
    input:
        bins=rules.move_genomes.output.dir,
        unbinned=rules.move_unbinned.output.fa,
    output:
        "genomes/alignments/all_contigs.fa",
    log:
        "logs/genomes/alignments/concat_genomes.log",
    params:
        ext="fa",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        """
        (cat {input.bins}/*{params.ext} {input.unbinned} > {output}) 1>{log} 2>&1
        """


# Skip indexing becuase it hard sets k-mer size and other params which we need to be flexible for 
# different input read types. I.e., optimal k-mer size varies between Nanopore vs. PacBio vs. Illumina
# reads. Thus having it hard set for all samples will lead to reduced accuray results for some samples.
rule align_reads_to_genomes:
    input:
        unpack(lambda wc: get_quality_controlled_reads(wc, as_dict=True)),
        target=rules.concat_genomes.output,
    output:
        "genomes/alignments/bams/{sample}.bam",
    params:
        command = lambda wc, input, output, threads, resources: align_reads_command(
            wc, input, output, threads, resources
        ),
    log:
        "logs/genomes/alignments/bams/{sample}_map.log",
    benchmark:
        "benchmarks/genomes/alignments/bams/{sample}_map.tsv",
    conda:
        "../envs/minimap.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "mapping", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule mapping_stats_genomes:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
    output:
        "genomes/alignments/stats/{sample}.stats",
    log:
        "logs/genomes/alignments/stats/{sample}_stats.log",
    threads: lambda wc: get_resource(wc, None, 1, "mapping_stats_genomes", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "account"),
    wrapper:
        "v1.19.0/bio/samtools/stats"


checkpoint multiqc_mapping_genome:
    input:
        expand("genomes/alignments/stats/{sample}.stats", sample=SAMPLES),
    output:
        "reports/quantify_genomes_mapping_results.html",
    log:
        "logs/genomes/alignment/multiqc.log",
    threads: lambda wc: get_resource(wc, None, 1, "multiqc_mapping_genome", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_genome", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_genome", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_genome", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_genome", "account"),
    wrapper:
        "v3.3.6/bio/multiqc"


rule mapping_coverm_coverage:
    input:
        g2c="genomes/clustering/{grouping}.genome2contig.tsv",
        bams=expand("genomes/alignments/bams/{sample}.bam", sample=SAMPLES),
    output:
        cov="genomes/coverage/{grouping}.coverage.tsv.gz",
        read_stats="genomes/coverage/{grouping}.read_stats.tsv",
    params:
        extra=config["coverm_params"],
        stats=config["coverm_stats"],
    log:
        general="logs/coverage/{grouping}.coverage.log",
        coverm="logs/coverage/{grouping}.coverm.log",
    threads: lambda wc: get_resource(wc, None, 1, "mapping_coverm_coverage", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "account"),
    conda:
        "../envs/coverm.yaml"
    shell:
        """
        (
        coverm genome \\
            {params.extra} \\
            --output-format sparse --methods {params.stats} \\
            --genome-definition {input.g2c} \\
            -t {threads} \\
            -b {input.bams} \\
          2>{log.coverm} \\
          | sed -e 's/\.coordSorted//g' \\
          | gzip -c \\
          1> {output.cov}; \\
        cat {log.coverm} \\
          | grep 'In sample' \\
          | sed -e "s/.* '\([^']*\).*found \([^ ]*\) reads mapped out of \([^ ]*\) total (\(.*\))/\\1\\t\\2\\t\\3\\t\\4/" \\
          | sort \\
          | awk 'BEGIN{{print "sample_id\\tmapped_reads\\ttotal_reads\\tpercent_mapped"}}{{print}}' \\
          1> {output.read_stats}
        ) 1>{log.general} 2>&1
        """


# Used in Snakefile
def get_coverm_files(wildcards):
    valid_files = []
    for grouping in ['prokaryotic', 'eukaryotic', 'viral', 'plasmid', 'unbinned', 'all', 'mags']:
        file_name = f"genomes/clustering/{grouping}.genome2contig.tsv"
        if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
            valid_files.append(f"genomes/coverage/{grouping}.coverage.tsv.gz")
    return valid_files



