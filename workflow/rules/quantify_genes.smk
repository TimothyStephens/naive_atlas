

#########################################
####                                 ####
#### Gene mapping and quantification ####
####                                 ####
#########################################

localrules:
    concat_cds,

rule concat_cds:
    input:
        binned_genes=lambda wc:   checkpoints.move_genome_predicted_genes.get(**wc).output.outdir,
        unbinned_genes=lambda wc: checkpoints.move_unbinned_predicted_genes.get(**wc).output.outdir,
    output:
        fa="genomes/alignments/all_cds.fa",
    log:
        "logs/genomes/genes/alignments/concat_cds.log",
    params:
        binned_fna_files=expand("genomes/genes/genomes/{genome}.fna",
            genome=get_all_output_predicted_genes('')),
        unbinned_fna_files=expand("genomes/genes/unbinned/{genome}.fna",
            genome=get_all_output_predicted_genes_unbinned('')),
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        """
        (
        cat {params.binned_fna_files} {params.unbinned_fna_files} > {output.fa}
        ) 1>{log} 2>&1
        """


# Skip indexing becuase it hard sets k-mer size and other params which we need to be flexible for 
# different input read types. I.e., optimal k-mer size varies between Nanopore vs. PacBio vs. Illumina
# reads. Thus having it hard set for all samples will lead to reduced accuray results for some samples.
rule align_reads_to_cds:
    input:
        unpack(lambda wc: get_quality_controlled_reads(wc, as_dict=True)),
        target=rules.concat_cds.output,
    output:
        temp("genomes/genes/alignments/bams/{sample}.bam"),
    params:
        command = lambda wc, input, output, threads, resources: align_reads_command(
            wc, input, output, threads, resources
        ),
    log:
        "logs/genomes/genes/alignments/bams/{sample}_map.log",
    benchmark:
        "benchmarks/genomes/genes/alignments/bams/{sample}_map.tsv",
    conda:
        "../envs/minimap.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "align_reads_to_cds", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "align_reads_to_cds", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "align_reads_to_cds", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "align_reads_to_cds", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "align_reads_to_cds", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """

rule mapping_stats_cds:
    input:
        bam="genomes/genes/alignments/bams/{sample}.bam",
    output:
        "genomes/genes/alignments/stats/{sample}.stats",
    log:
        "logs/genomes/genes/alignments/stats/{sample}_stats.log",
    threads: lambda wc: get_resource(wc, None, 1, "mapping_stats_cds", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_cds", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_cds", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_cds", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_cds", "account"),
    wrapper:
        "v1.19.0/bio/samtools/stats"


checkpoint multiqc_mapping_cds:
    input:
        expand("genomes/genes/alignments/stats/{sample}.stats", sample=SAMPLES),
    output:
        "reports/quantify_cds_mapping_results.html",
    log:
        "logs/genomes/genes/alignment/multiqc.log",
    threads: lambda wc: get_resource(wc, None, 1, "multiqc_mapping_cds", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_cds", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_cds", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_cds", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "multiqc_mapping_cds", "account"),
    wrapper:
        "v3.3.6/bio/multiqc"


rule mapping_coverm_coverage_cds:
    input:
        bams=expand("genomes/genes/alignments/bams/{sample}.bam", sample=SAMPLES),
    output:
        cov="genomes/genes/coverage/all_cds.coverage.tsv.gz",
        read_stats="genomes/genes/coverage/all_cds.read_stats.tsv",
    params:
        extra=config["coverm_params"],
        stats=config["coverm_stats"],
    log:
        general="logs/genomes/genes/coverage/all_cds.coverage.log",
        coverm="logs/genomes/genes/coverage/all_cds.coverm.log",
    benchmark:
        "benchmarks/genomes/genes/coverage/all_cds.coverm.tsv"
    threads: lambda wc: get_resource(wc, None, 1, "mapping_coverm_coverage_cds", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage_cds", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage_cds", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage_cds", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage_cds", "account"),
    conda:
        "../envs/coverm.yaml"
    shell:
        """
        (
        coverm contig \\
            {params.extra} \\
            --output-format sparse --methods {params.stats} \\
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


