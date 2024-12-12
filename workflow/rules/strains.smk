


rule instrain_profile:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
        genomes=rules.concat_genomes.output,
        # genes=lambda wc: get_all_genes(wc, extension=".fna"),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("Intermediate/strains/{sample}"),
    threads: config["simplejob_threads"]
    params:
        extra=config.get("instrain_profile_extra", ""),
    log:
        "logs/genomes/strains/profile/{sample}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/genomes/strains/profile/{sample}.tsv"
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        #" cat {input.genes} > {resources.tmpdir}/all_genome_genes.fna 2> {log} "
        #" ; "
        "inStrain profile "
        " {input.bam} {input.genomes} "
        " -o {output} "
        " -p {threads} "

        " -s {input.scaffold_to_genome} "
        " --database_mode "
        " {params.extra} &>> {log}"
        #" -g {resources.tmpdir}/all_genome_genes.fna "


rule instrain_compare:
    input:
        profiles=expand("Intermediate/strains/{sample}", sample=SAMPLES),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("genomes/strains/comparison"),
    threads: config["simplejob_threads"]
    params:
        extra=config.get("instrain_compare_extra", ""),
    log:
        "logs/genomes/strains/compare.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/genomes/strains/compare.tsv"
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        "inStrain compare "
        " --input {input.profiles} "
        " -o {output} "
        " -p {threads} "
        " -s {input.scaffold_to_genome} "
        " --database_mode "
        " {params.extra} &> {log}"


# usage: inStrain compare -i [INPUT [INPUT ...]] [-o OUTPUT] [-p PROCESSES] [-d]
#                         [-h] [--version] [-s [STB [STB ...]]] [-c MIN_COV]
#                         [-f MIN_FREQ] [-fdr FDR] [--database_mode]
#                         [--breadth BREADTH] [-sc SCAFFOLDS] [--genome GENOME]
#                         [--store_coverage_overlap]
#                         [--store_mismatch_locations]
#                         [--include_self_comparisons] [--skip_plot_generation]
#                         [--group_length GROUP_LENGTH] [--force_compress]
#                         [-ani ANI_THRESHOLD] [-cov COVERAGE_TRESHOLD]
#                         [--clusterAlg {ward,single,complete,average,weighted,median,centroid}]
