import os


#######################
####               ####
####    EGG NOG    ####
####               ####
#######################

# output with wildcards "{folder}/{prefix}.emapper.tsv"

rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa="genomes/genomes/{genome}.faa",
    output:
        temp(
            "Intermediate/genecatalog/annotation/eggNOG/{genome}.emapper.seed_orthologs"
        ),
    params:
        data_dir=EGGNOG_DIR,
        prefix=lambda wc, output: output[0].replace(".emapper.seed_orthologs", ""),
    resources:
        mem=config["simplejob_memory"],
    threads: config["simplejob_threads"]
    shadow:
        "minimal"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/genecatalog/annotation/eggnog/{genome}_homology_search_diamond.log",
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> {log}
        """


def calculate_mem_eggnog():
    return 2 * config["simplejob_memory"] + (
        37 if config["eggNOG_use_virtual_disk"] else 0
    )


rule eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed=rules.eggNOG_homology_search.output,
    output:
        temp("Intermediate/genecatalog/annotation/eggNOG/{genome}.emapper.annotations"),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else EGGNOG_DIR
        ),
        prefix=lambda wc, output: output[0].replace(".emapper.annotations", ""),
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: config["simplejob_threads"]
    resources:
        mem=calculate_mem_eggnog(),
    shadow:
        "minimal"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/genecatalog/annotation/eggnog/{genome}_annotate_hits_table.log",
    shell:
        """
        if [ {params.copyto_shm} == "t" ] ; then
            # Check if the files exist before copying
            if [ ! -e "{params.data_dir}/eggnog.db" ]; then
                cp {EGGNOG_DIR}/eggnog.db {params.data_dir}/eggnog.db 2> {log}
            else
                echo "File {params.data_dir}/eggnog.db already exists. Skipping copy." >> {log}
            fi

            if [ ! -e "{params.data_dir}/eggnog_proteins.dmnd" ]; then
                cp {EGGNOG_DIR}/eggnog_proteins.dmnd {params.data_dir}/eggnog_proteins.dmnd 2>> {log}
            else
                echo "File {params.data_dir}/eggnog_proteins.dmnd already exists. Skipping copy." >> {log}
            fi
        fi

        emapper.py --annotate_hits_table {input.seed} --no_file_comments \
          --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2>> {log}

        """


def get_all_gene_eggnog(wildcards):
    all_genomes = get_all_genomes(wildcards)
    return expand(rules.eggNOG_annotation.output, genome=all_genomes)

rule combine_egg_nog_annotations:
    input:
        get_all_gene_eggnog,
    output:
        parquet="genomes/annotations/genes/eggNOG.parquet",
        tsv="genomes/annotations/genes/eggNOG.tsv.gz",
    log:
        "logs/genomes/annotations/genes/eggNOG/combine.log",
    resources:
        time=config["simplejob_runtime"],
    run:
        try:
            import pandas as pd

            Tables = [
                pd.read_csv(file, index_col=None, header=0, sep="\t")
                for file in input
            ]

            combined = pd.concat(Tables, axis=0)

            del Tables

            combined.columns = EGGNOG_HEADER
            #combined["Seed_evalue"] = combined["Seed_evalue"].astype("bytes")
            #combined["Seed_Score"] = combined["Seed_Score"].astype("bytes")

            combined.to_parquet(output["parquet"], index=False)
            combined.to_csv(output["tsv"], sep='\t', index=False)
        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e





######################
####              ####
####     DRAM     ####
####              ####
######################

rule DRAM_annotation:
    input:
        faa="genomes/genomes/{genome}.faa",
        config=get_dram_config,
    output:
        annotations=temp(
            "Intermediate/genecatalog/annotation/dram/{genome}/annotations.tsv"
        ),
        genes=temp("Intermediate/genecatalog/annotation/dram/{genome}/genes.faa"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    params:
        extra=config.get("dram_extra", ""),
        outdir=lambda wc, output: Path(output[0]).parent,
    log:
        "logs/Genecatalog/annotation/dram/{genome}.log",
        "logs/Genecatalog/annotation/dram/{genome}.logfile",
    shell:
        " rm -rf {params.outdir} &> {log[0]};"
        "\n"
        " DRAM.py annotate_genes "
        " --input_faa {input.faa}"
        " --config_loc {input.config} "
        " --output_dir {params.outdir} "
        " --threads {threads} "
        " {params.extra} "
        " --log_file_path {log[1]} "
        " --verbose &>> {log[0]}"


def get_all_gene_dram(wildcards):
    all_genomes = get_all_genomes(wildcards)
    return expand(rules.DRAM_annotation.output.annotations, genome=all_genomes)

rule combine_dram_genecatalog_annotations:
    input:
        get_all_gene_dram,
    output:
        directory("genomes/annotations/genes/dram"),
    resources:
        time=config["simplejob_runtime"],
    log:
        "logs/genomes/annotations/genes/dram/combine.log",
    script:
        "../scripts/combine_dram_gene_annotations.py"



