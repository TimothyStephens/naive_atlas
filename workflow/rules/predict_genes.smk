
localrules:
    copy_eukaryotic_genomes,

rule copy_eukaryotic_genomes:
    input:
        "genomes/genomes",
    output:
        directory("tmp/genes/eukaryotic"),
    log:
        "logs/gene_prediction/copy_eukaryotic_genomes.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources: 
        lambda wc, input, attempt: get_all_resources(wc, input, attempt, "localrule")
    shell:
        "mkdir -p {output} && cp {input}/MAG_eukaryotic* {output} > {log} 2>&1"


rule predict_eukaryotic_genes:
    input:
        genomes=rules.copy_eukaryotic_genomes.output,
    output:
        "samples/{sample}/binning/coverage/{sample_reads}.metabat_depth.txt",
    benchmark:
        "benchmarks/samples/{sample}/binning/coverage/{sample_reads}.txt"
    log:
        "logs/samples/{sample}/binning/coverage/{sample_reads}.log",
    conda:
        "../envs/metabat.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "predict_genes", "threads")
    resources: 
        lambda wc, input, attempt: get_all_resources(wc, input, attempt, "predict_genes")
    params:
        minid=config["cobinning_readmapping_id"] * 100,
    priority: 100
    shell:
        "jgi_summarize_bam_contig_depths "
        " --percentIdentity {params.minid} "
        " --outputDepth {output} "
        " {input} &> {log} "
