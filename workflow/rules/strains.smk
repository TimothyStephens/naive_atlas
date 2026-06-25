import glob

rule combine_and_convert_genes:
    input:
        binned_genes=lambda wc:   checkpoints.move_genome_predicted_genes.get(**wc).output.outdir,
        unbinned_genes=lambda wc: checkpoints.move_unbinned_predicted_genes.get(**wc).output.outdir,
    params:
        gff=glob.glob("genome/genes/*/*.gff3"),
        fna=glob.glob("genome/genes/*/*.fna"),
    output:
        fna="genomes/strains/combined_genes.fna",
        skipped="genomes/strains/skipped_multiexon_genes.txt",
        mapping="genomes/strains/gene_id_mapping.tsv",
    log:
        "logs/genomes/strains/combine_and_convert_genes.log",
    benchmark:
        "benchmarks/genomes/strains/combine_and_convert_genes.tsv",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    script:
        "../scripts/gff_to_prodigal.py"


rule instrain_profile:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
        genomes=rules.concat_genomes.output,
        scaffold_to_genome="genomes/clustering/mags.contig2genome.tsv",
        genes=rules.combine_and_convert_genes.output.fna,
    output:
        directory("genomes/strains/profiles/{sample}"),
    threads: lambda wc: get_resource(wc, None, 1, "instrain_profile", "threads")
    params:
        extra=config["instrain_profile_extra"],
    log:
        "logs/genomes/strains/profiles/{sample}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "benchmarks/genomes/strains/profiles/{sample}.tsv",
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
            -g {input.genes} \\
            -o {output} \\
            -p {threads} \\
            -s {input.scaffold_to_genome} \\
            --database_mode {params.extra}
        ) 1>{log} 2>&1
        """


rule instrain_compare:
    input:
        profiles=expand("genomes/strains/profiles/{sample}", sample=SAMPLES),
        scaffold_to_genome="genomes/clustering/mags.contig2genome.tsv",
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
            --database_mode {params.extra}
        ) 1>{log} 2>&1
        """


