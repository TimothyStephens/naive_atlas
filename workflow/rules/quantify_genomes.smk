




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
        exit(1)

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
        exit(1)

    return genomes


rule get_contig2genomes:
    input:
        genome_dir,
    output:
        c2g="genomes/clustering/contig2genome.tsv",
        g2c="genomes/clustering/genome2contig.tsv",
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
    shell:
        "cat {input}/*{params.ext} > {output}"


if config["genome_aligner"] == "minimap":

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
        threads: 3
        resources:
            mem=config["simplejob_memory"],
        wrapper:
            "v1.19.0/bio/minimap2/index"

    rule align_reads_to_genomes:
        input:
            target=rules.index_genomes.output,
            query=get_quality_controlled_reads,
        output:
            "genomes/alignments/bams/{sample}.bam",
        log:
            "logs/genomes/alignments/{sample}_map.log",
        params:
            extra="-x sr",
            sorting="coordinate",
        threads: config["simplejob_threads"]
        resources:
            mem=config["simplejob_memory"],
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

elif config["genome_aligner"] == "bwa":

    rule index_genomes:
        input:
            rules.concat_genomes.output,
            timestamp=genome_dir,
        output:
            multiext("ref/genomes", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        log:
            "logs/genomes/alignments/bwa_index.log",
        threads: 4
        resources:
            mem=config["simplejob_memory"],
        wrapper:
            "v1.19.0/bio/bwa-mem2/index"

    rule align_reads_to_genomes:
        input:
            idx=rules.index_genomes.output,
            reads=get_quality_controlled_reads,
        output:
            "genomes/alignments/bams/{sample}.bam",
        log:
            "logs/genomes/alignments/{sample}_bwa.log",
        params:
            extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
            sort="samtools",
            sort_order="coordinate",
        threads: config["simplejob_threads"]
        resources:
            mem=config["simplejob_memory"],
            mem_mb=config["simplejob_memory"] * 1000,
        wrapper:
            "v1.19.0/bio/bwa-mem2/mem"

else:
    raise Exception(
        "'genome_aligner' not understood, it should be 'minimap' or 'bwa', not '{genome_aligner}'. check config file".format(
            **config
        )
    )


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
    shell:
        "mv {input} {output} > {log}"


rule mapping_stats_genomes:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
    output:
        "genomes/alignments/stats/{sample}.stats",
    log:
        "logs/genomes/alignments/{sample}_stats.log",
    resources:
        mem=config["simplejob_memory"],
    wrapper:
        "v1.19.0/bio/samtools/stats"


rule multiqc_mapping_genome:
    input:
        expand("genomes/alignments/stats/{sample}.stats", sample=SAMPLES),
    output:
        "reports/genome_mapping/results.html",
    log:
        "logs/genomes/alignment/multiqc.log",
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
    threads: config["simplejob_threads"]
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



