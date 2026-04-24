

#####################################
####                             ####
#### Genome mapping and coverage ####
####                             ####
#####################################

def get_genome_dir():
    if ("genome_dir" in config) and (config["genome_dir"] is not None):
        genome_dir = config["genome_dir"]
        assert os.path.exists(genome_dir), f"{genome_dir} Doesn't exists"
        
        logger.info(f"Set genomes from {genome_dir}.")
        
        # check if genomes are present
        genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fa")).genome
        
        if len(genomes) == 0:
            logger.error(f"No genomes found with fa extension in {genome_dir} ")
            exit(1)
    
    else:
        genome_dir = "genomes/genomes"
    
    return genome_dir


genome_dir = get_genome_dir()


def get_all_genomes(wildcards):
    global genome_dir
    # check if genomes are present
    genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fa")).genome
    
    if len(genomes) == 0:
        logger.error(
            f"No genomes found with fa extension in {genome_dir} "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        #exit(1)

    return genomes


def get_all_unbinned(wildcards):
    # check if genomes are present
    genomes = glob_wildcards(os.path.join("genomes/unbinned", "{genome}.fa")).genome

    if len(genomes) == 0:
        logger.error(
            f"No genomes found with fasta extension in genomes/genomes/unbinned "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        #exit(1)

    return genomes


localrules:
    get_contig2genomes,

rule get_contig2genomes:
    input:
        genome_dir,
    output:
        c2g="genomes/clustering/contig2genome.tsv",
        g2c="genomes/clustering/genome2contig.tsv",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        from glob import glob

        fasta_files = glob(input[0] + "/*.fa")

        with open(output["c2g"], "w") as out_c2g, open(output["g2c"], "w") as out_g2c:
            for fasta in fasta_files:
                bin_name, ext = os.path.splitext(os.path.split(fasta)[-1])
                # if gz remove also fasta extension
                if ext == ".gz":
                    bin_name = os.path.splitext(bin_name)[0]

                    # write names of contigs in mapping file
                with open(fasta) as f:
                    for line in f:
                        if line[0] == ">":
                            header = line[1:].strip().split()[0]
                            out_c2g.write(f"{header}\t{bin_name}\n")
                            out_g2c.write(f"{bin_name}\t{header}\n")

# alternative way to get to contigs2genomes for quantification with external genomes
ruleorder: get_contig2genomes > rename_genomes


### Quantification

localrules:
    concat_genomes,

rule concat_genomes:
    input:
        genome_dir,
    output:
        "genomes/alignments/all_contigs.fa",
    params:
        ext="fa",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        "cat {input}/*{params.ext} > {output}"


rule index_genomes:
    input:
        target=rules.concat_genomes.output,
        timestamp=genome_dir,
    output:
        "ref/genomes.mmi",
    log:
        "logs/genomes/alignmentsindex.log",
    params:
        index_size="12G",
    threads: lambda wc: get_resource(wc, None, 1, "mapping", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping", "account"),
    wrapper:
        "v3.13.4/bio/minimap2/index"


rule align_reads_to_genomes:
    input:
        unpack(lambda wc: get_quality_controlled_reads(wc, as_dict=True)),
        target=rules.index_genomes.output,
    output:
        "genomes/alignments/bams/{sample}.bam",
    params:
        command = lambda wc, input, output, threads, resources: align_reads_command(
            wc, input, output, threads, resources
        ),
    log:
        "logs/genomes/alignments/{sample}_map.log",
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
        ({params.command}) > {log} 2>&1
        """



# path change for bam file
localrules:
    move_old_bam,


ruleorder: move_old_bam > align_reads_to_genomes


rule move_old_bam:
    input:
        "genomes/alignments/{sample}.bam",
    output:
        "genomes/alignments/bams/{sample}.bam",
    log:
        "logs/genomes/alignments/{sample}_move.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        "mv {input} {output} > {log}"


rule mapping_stats_genomes:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
    output:
        "genomes/alignments/stats/{sample}.stats",
    log:
        "logs/genomes/alignments/{sample}_stats.log",
    threads: lambda wc: get_resource(wc, None, 1, "mapping_stats_genomes", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_stats_genomes", "account"),
    wrapper:
        "v1.19.0/bio/samtools/stats"


rule multiqc_mapping_genome:
    input:
        expand("genomes/alignments/stats/{sample}.stats", sample=SAMPLES),
    output:
        "reports/genome_mapping/results.html",
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
        g2c="genomes/clustering/genome2contig.tsv",
        bams=expand("genomes/alignments/bams/{sample}.bam", sample=SAMPLES),
    output:
        cov="genomes/coverage/coverage.tsv.gz",
        read_stats="genomes/coverage/read_stats.tsv",
    params:
        extra=config["coverm_params"],
        stats=config["coverm_stats"],
    log:
        general="logs/coverage/coverage.log",
        coverm="logs/coverage/coverm.log",
    threads: lambda wc: get_resource(wc, None, 1, "mapping_coverm_coverage", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "mapping_coverm_coverage", "account"),
    conda:
        "../envs/coverm.yaml"
    shell:
        "("
        "coverm genome"
        "  {params.extra}"
        "  --output-format sparse --methods {params.stats}"
        "  --genome-definition {input.g2c}"
        "  -t {threads}"
        "  -b {input.bams}"
        " 2>{log.coverm}"
        " | sed -e 's/\.coordSorted//g'"
        " | gzip -c"
        " > {output.cov}; "
        "cat {log.coverm}"
        " | grep 'In sample'"
        " | sed -e \"s/.* '\([^']*\).*found \([^ ]*\) reads mapped out of \([^ ]*\) total (\(.*\))/\\1\\t\\2\\t\\3\\t\\4/\""
        " | sort"
        " | awk 'BEGIN{{print \"sample_id\\tmapped_reads\\ttotal_reads\\tpercent_mapped\"}}{{print}}'"
        " > {output.read_stats}"
        ")"
        " 1>{log.general} 2>&1"


