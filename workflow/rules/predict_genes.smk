


rule copy_eukaryotic_genomes:
    input:
        GENOME_DIR,
    output:
        directory("tmp/genes/eukaryotic"),
    log:
        "logs/gene_prediction/copy_eukaryotic_genomes.log",
    shell:
        "mkdir -p {output} && cp {input}/MAG_eukaryotic* {output}"


rule predict_eukaryotic_genes:
    input:
        genomes=rules.copy_eukaryotic_genomes.output,
    output:
        "{sample}/binning/coverage/{sample_reads}.metabat_depth.txt",
    benchmark:
        "{sample}/logs/benchmarks/binning/coverage/{sample_reads}.txt"
    log:
        "{sample}/logs/binning/coverage/{sample_reads}.log",
    conda:
        "../envs/metabat.yaml"
    threads: config["simplejob_threads"]  # multithreaded trough OMP_NUM_THREADS
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    params:
        minid=config["cobinning_readmapping_id"] * 100,
    priority: 100
    shell:
        "jgi_summarize_bam_contig_depths "
        " --percentIdentity {params.minid} "
        " --outputDepth {output} "
        " {input} &> {log} "



