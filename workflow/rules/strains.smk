
rule instrain_profile:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
        genomes=rules.concat_genomes.output,
        # genes=lambda wc: get_all_genes(wc, extension=".fna"),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("Intermediate/strains/{sample}"),
    threads: lambda wc: get_resource(wc, None, 1, "instrain_profile", "threads")
    params:
        extra=config["instrain_profile_extra"],
    log:
        "logs/genomes/strains/profile/{sample}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "benchmarks/genomes/strains/profile/{sample}.tsv",
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_profile", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_profile", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_profile", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_profile", "account"),
    shell:
        """
        (
        inStrain profile \\
            {input.bam} {input.genomes} \\
            -o {output} \\
            -p {threads} \\
            -s {input.scaffold_to_genome} \\
            --database_mode \\
            {params.extra}
        ) 1>{log} 2>&1
        """


rule instrain_compare:
    input:
        profiles=expand("Intermediate/strains/{sample}", sample=SAMPLES),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("genomes/strains/comparison"),
    threads: lambda wc: get_resource(wc, None, 1, "instrain_compare", "threads")
    params:
        extra=config["instrain_compare_extra"],
    log:
        "logs/genomes/strains/compare.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "benchmarks/genomes/strains/compare.tsv",
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_compare", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_compare", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_compare", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "instrain_compare", "account"),
    shell:
        """
        (
        inStrain compare \\
            --input {input.profiles} \\
            -o {output} \\
            -p {threads} \\
            -s {input.scaffold_to_genome} \\
            --database_mode \\
            {params.extra}
        ) 1>{log} 2>&1
        """


