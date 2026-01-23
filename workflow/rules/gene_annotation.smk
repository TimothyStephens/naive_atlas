import os


#######################
####               ####
####    EGG NOG    ####
####               ####
#######################

# output with wildcards "{folder}/{prefix}.emapper.tsv"

rule gene_eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa="genomes/genes/{dataset}/{genome}.faa",
    output:
        temp(
            "Intermediate/genecatalog/annotations/{dataset}/genes/eggNOG/{genome}.emapper.seed_orthologs"
        ),
    params:
        data_dir=EGGNOG_DIR,
        prefix=lambda wc, output: output[0].replace(".emapper.seed_orthologs", ""),
    resources:
        mem=config["simplejob_memory"],
    threads: config["simplejob_threads"]
    #conda:
    #    "../envs/eggNOG.yaml"
    container:
        "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1"
    log:
        "logs/genecatalog/annotations/{dataset}/genes/eggnog/{genome}_homology_search_diamond.log",
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override &> {log}
        """


def calculate_mem_eggnog():
    return 2 * config["simplejob_memory"] + (
        37 if config["eggNOG_use_virtual_disk"] else 0
    )


rule gene_eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed=rules.gene_eggNOG_homology_search.output,
    output:
        temp("Intermediate/genecatalog/annotations/{dataset}/genes/eggNOG/{genome}.emapper.annotations"),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else EGGNOG_DIR
        ),
        prefix=lambda wc, output: output[0].replace(".emapper.annotations", ""),
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: config["simplejob_threads"]
    resources:
        mem=calculate_mem_eggnog(),
    #conda:
    #    "../envs/eggNOG.yaml"
    container:
        "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1"
    log:
        "logs/genecatalog/annotations/{dataset}/genes/eggnog/{genome}_annotate_hits_table.log",
    shell:
        """
        if [ {params.copyto_shm} == "t" ] ; then
            # Check if the files exist before copying
            if [ ! -e "{params.data_dir}/eggnog.db" ]; then
                cp {EGGNOG_DIR}/eggnog.db {params.data_dir}/eggnog.db &> {log}
            else
                echo "File {params.data_dir}/eggnog.db already exists. Skipping copy." &>> {log}
            fi

            if [ ! -e "{params.data_dir}/eggnog_proteins.dmnd" ]; then
                cp {EGGNOG_DIR}/eggnog_proteins.dmnd {params.data_dir}/eggnog_proteins.dmnd &>> {log}
            else
                echo "File {params.data_dir}/eggnog_proteins.dmnd already exists. Skipping copy." &>> {log}
            fi
        fi

        emapper.py --annotate_hits_table {input.seed} --no_file_comments \
          --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} &>> {log}

        """


def get_all_gene_eggnog(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return expand(rules.gene_eggNOG_annotation.output,
            dataset=wildcards.dataset, genome=all_genomes)

rule combine_gene_egg_nog_annotations:
    input:
        get_all_gene_eggnog,
    output:
        parquet="genomes/annotations/{dataset}/genes/eggNOG.parquet",
        tsv="genomes/annotations/{dataset}/genes/eggNOG.tsv.gz",
    log:
        "logs/genomes/annotations/{dataset}/genes/eggNOG/combine.log",
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
            combined["Seed_evalue"] = combined["Seed_evalue"].astype("bytes")
            combined["Seed_Score"] = combined["Seed_Score"].astype("bytes")
            combined = combined.astype(str)
            
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

rule gene_DRAM_annotation:
    input:
        faa="genomes/genes/{dataset}/{genome}.faa",
        config=get_dram_config,
    output:
        annotations=temp(
            "Intermediate/genecatalog/annotations/{dataset}/genes/dram/{genome}/annotations.tsv"
        ),
        genes=temp("Intermediate/genecatalog/annotations/{dataset}/genes/dram/{genome}/genes.faa"),
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
        "logs/Genecatalog/annotations/{dataset}/genes/dram/{genome}.log",
        "logs/Genecatalog/annotations/{dataset}/genes/dram/{genome}.logfile",
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
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return expand(rules.gene_DRAM_annotation.output.annotations,
                dataset=wildcards.dataset, genome=all_genomes)

rule combine_gene_dram_genecatalog_annotations:
    input:
        get_all_gene_dram,
    output:
        directory("genomes/annotations/{dataset}/genes/dram"),
    resources:
        time=config["simplejob_runtime"],
    log:
        "logs/genomes/annotations/{dataset}/genes/dram/combine.log",
    script:
        "../scripts/combine_dram_gene_annotations.py"





######################
####              ####
####   MMSEQS2    ####
####              ####
######################

rule gene_mmseqs2_annotation:
    input:
        faa="genomes/genes/{dataset}/{genome}.faa",
        database=rules.mmseqs2_download.output.database,
    output:
        results="genomes/annotations/{dataset}/genes/mmseqs2/{genome}.faa.mmseqs2_{database_name}.m4.gz",
        tmp=temp(directory("Intermediate/annotations/{dataset}/genes/mmseqs2/{genome}.faa.mmseqs2_{database_name}.tmp")),
    params:
        mmseqs2_opts=config["mmseqs2_opts"],
        mem=int(config["mmseqs2_memory"]*0.8),
        results="genomes/annotations/{dataset}/genes/mmseqs2/{genome}.faa.mmseqs2_{database_name}.m4",
    threads: config["simplejob_threads"]
    resources:
        mem=config["mmseqs2_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "logs/genomes/annotations/{dataset}/genes/mmseqs2/{database_name}/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/genes/mmseqs2/{database_name}/{genome}.tsv"
    shell:
        """
        (
        mmseqs easy-search \
            --threads {threads} \
            --split-memory-limit {params.mem}G \
            --compressed 1 \
            --format-mode 4 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,taxid,taxname,taxlineage,theader \
            {params.mmseqs2_opts} \
            {input.faa} \
            {input.database} \
            {params.results} \
            {output.tmp} \
          && pigz -11 -p {threads} {params.results}
        ) &> {log}
        """


def get_all_gene_mmseqs2_annotation(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return(expand(rules.gene_mmseqs2_annotation.output.results,
                        dataset=wildcards.dataset,
                        database_name=config["mmseqs2_database_name"],
                        genome=all_genomes
            )
    )

localrules:
    all_mmseqs2,

rule all_mmseqs2:
    input:
        get_all_gene_mmseqs2_annotation,
    output:
        touch("genomes/annotations/{dataset}/genes/mmseqs2/finished"),



