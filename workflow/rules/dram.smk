DBDIR = config["database_dir"]


def get_dram_config(wildcards):
    old_dram_path = f"{DBDIR}/Dram"
    if Path(old_dram_path).exists():
        logger.error(
            f"Detected an old database for DRAM in {old_dram_path}. You can delete it."
        )

    return config.get("dram_config_file", f"{DBDIR}/DRAM/DRAM.config")


rule DRAM_annotate:
    input:
        fasta="genomes/{dataset}/{genome}.fa",
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/{dataset}/dram/intermediate_files/{genome}"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    params:
        extra=config.get("dram_extra", ""),
        min_contig_size=config.get("minimum_contig_length", "1000"),
    log:
        "logs/annotations/{dataset}/dram/run_dram/{genome}.log",
    benchmark:
        "logs/benchmarks/annotations/{dataset}/dram/run_dram/{genome}.tsv"
    shell:
        " DRAM.py annotate "
        " --config_loc {input.config} "
        " --input_fasta {input.fasta}"
        " --output_dir {output.outdir} "
        " --threads {threads} "
        " --min_contig_size {params.min_contig_size} "
        " {params.extra} "
        " --verbose &> {log}"
        #" --gtdb_taxonomy {input.gtdb_dir}/{params.gtdb_file} "
        #" --checkm_quality {input.checkm} "


def get_all_dram(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return expand(rules.DRAM_annotate.output.outdir,
            dataset=wildcards.dataset, genome=all_genomes)


localrules:
    concat_annotations,

rule concat_annotations:
    input:
        get_all_dram,
    output:
        "genomes/annotations/{dataset}/dram/annotations.tsv",
    resources:
        time=config["simplejob_runtime"],
    run:
        from utils import io

        for i, annotation_file in enumerate(["annotations.tsv"]):
            input_files = [
                os.path.join(dram_folder, annotation_file) for dram_folder in input
            ]

            io.pandas_concat(
                input_files, output[i], sep="\t", index_col=0, axis=0, disk_based=True
            )


rule DRAM_destill:
    input:
        rules.concat_annotations.output,
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/{dataset}/dram/distil"),
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/annotations/{dataset}/dram/distil.log",
    shell:
        " DRAM.py distill "
        " --config_loc {input.config} "
        " --input_file {input[0]}"
        " --output_dir {output} "
        "  &> {log}"


rule get_all_modules:
    input:
        annotations="genomes/annotations/{dataset}/dram/annotations.tsv",
        config=get_dram_config,
    output:
        "genomes/annotations/{dataset}/dram/kegg_modules.tsv",
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/annotations/{dataset}/dram/get_all_modules.log",
    script:
        "../scripts/DRAM_get_all_modules.py"


rule dram:
    input:
        "genomes/annotations/{dataset}/dram/distil",
        "genomes/annotations/{dataset}/dram/kegg_modules.tsv",
    output:
        touch("genomes/annotations/{dataset}/dram/finished"),


