import os


#######################
####               ####
####    EGG NOG    ####
####               ####
#######################

# output with wildcards "{folder}/{prefix}.emapper.tsv"

rule gene_eggNOG_homology_search:
    input:
        eggnog_db_files=rules.download_eggNOG_files.output.files,
        faa="genomes/genes/{dataset}/{genome}.faa",
    output:
        seed=temp(
            "Intermediate/genecatalog/annotations/{dataset}/genes/eggNOG/{genome}.emapper.seed_orthologs"
        ),
        hits=temp(
            "Intermediate/genecatalog/annotations/{dataset}/genes/eggNOG/{genome}.emapper.hits"
        ),
    params:
        data_dir=rules.download_eggNOG_files.output.dir,
        prefix=lambda wc, output: output[0].replace(".emapper.seed_orthologs", ""),
    threads: lambda wc: get_resource(wc, None, 1, "gene_annot_eggnog", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "account"),
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


rule gene_eggNOG_annotation:
    input:
        eggnog_db_files=rules.download_eggNOG_files.output.files,
        seed=rules.gene_eggNOG_homology_search.output.seed,
    output:
        temp("Intermediate/genecatalog/annotations/{dataset}/genes/eggNOG/{genome}.emapper.annotations"),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else rules.download_eggNOG_files.output.dir
        ),
        prefix=lambda wc, output: output[0].replace(".emapper.annotations", ""),
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: lambda wc: get_resource(wc, None, 1, "gene_annot_eggnog", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "account"),
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
    """
    Make sure we only request annotation files for *.faa files with sequences in them.
    """
    # Make sure we have finished moving the final gene prediction files.
    checkpoint_output = checkpoints.move_genome_predicted_genes.get(**wildcards).output.outdir
    
    # Look for all *.faa files in the `genomes/genes` directory
    # Expect: genomes/genes/{dataset}/{genome}.faa
    FAA_FILES = glob_wildcards(f"genomes/genes/{wildcards.dataset}/{{genome}}.faa")
    
    valid_paths = []
    for genome in sorted(FAA_FILES.genome):
        path = f"genomes/genes/{wildcards.dataset}/{genome}.faa"
        
        if os.path.exists(path) and os.path.getsize(path) > 0:
            valid_paths.append(expand(rules.gene_eggNOG_annotation.output, dataset=wildcards.dataset, genome=genome)[0])
        else:
            if config.get("debug", False):
                print(f"[DEBUG] Skipping eggNOG for {path}: File empty or missing.")
    
    return valid_paths


rule combine_gene_egg_nog_annotations:
    input:
        get_all_gene_eggnog,
    output:
        parquet="genomes/annotations/{dataset}/genes/eggNOG.parquet",
        tsv="genomes/annotations/{dataset}/genes/eggNOG.tsv.gz",
    log:
        "logs/genomes/annotations/{dataset}/genes/eggNOG/combine.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
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
####   MMSEQS2    ####
####              ####
######################

rule gene_mmseqs2_annotation:
    input:
        faa="genomes/genes/{dataset}/{genome}.faa",
        database=rules.mmseqs2_download.output.database,
    output:
        results="genomes/annotations/{dataset}/genes/mmseqs2_easy_search/{genome}.faa.mmseqs2_{database_name}.m4.gz",
        tmp=temp(directory("Intermediate/annotations/{dataset}/genes/mmseqs2_easy_search/{genome}.faa.mmseqs2_{database_name}.tmp")),
    params:
        mmseqs2_opts=config["mmseqs2_opts"],
        results="genomes/annotations/{dataset}/genes/mmseqs2_easy_search/{genome}.faa.mmseqs2_{database_name}.m4",
    threads: lambda wc: get_resource(wc, None, 1, "gene_annot_mmseqs2_easy_search", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_mmseqs2_easy_search", "mem_mb"),
        mem_gb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_mmseqs2_easy_search", "mem_gb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_mmseqs2_easy_search", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_mmseqs2_easy_search", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_mmseqs2_easy_search", "account"),
    container:
        "docker://ghcr.io/soedinglab/mmseqs2:18-8cc5c"
    log:
        "logs/genomes/annotations/{dataset}/genes/mmseqs2_easy_search/{database_name}/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/genes/mmseqs2_easy_search/{database_name}/{genome}.tsv"
    shell:
        """
        (
        /usr/local/bin/entrypoint easy-search \
            --threads {threads} \
            --split-memory-limit {resources.mem_gb}G \
            --format-mode 4 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,taxid,taxname,taxlineage,theader \
            {params.mmseqs2_opts} \
            {input.faa} \
            {input.database} \
            {params.results} \
            {output.tmp} \
          && gzip -9 {params.results}
        ) &> {log}
        """


def get_all_gene_mmseqs2_annotation(wildcards):
    """
    Make sure we only request annotation files for *.faa files with sequences in them.
    """
    # Make sure we have finished moving the final gene prediction files.
    checkpoint_output = checkpoints.move_genome_predicted_genes.get(**wildcards).output.outdir
    
    # Look for all *.faa files in the `genomes/genes` directory
    # Use an f-string to resolve 'wildcards.dataset' into 'unbinned' or 'genomes'
    search_pattern = f"genomes/genes/{wildcards.dataset}/{{genome}}.faa"
    
    # Now glob_wildcards only sees '{genome}'
    FAA_FILES = glob_wildcards(search_pattern)
    
    valid_paths = []
    for genome in FAA_FILES.genome:
        path = f"genomes/genes/{wildcards.dataset}/{genome}.faa"
        
        # Check if file exists and has content
        if os.path.exists(path) and os.path.getsize(path) > 0:
            # Append the expected output path for the annotation rule
            valid_paths.append(
                expand(rules.gene_mmseqs2_annotation.output.results,
                       dataset=wildcards.dataset,
                       database_name=config["mmseqs2_database_name"],
                       genome=genome)[0]
            )
        else:
            if config.get("debug", False):
                print(f"[DEBUG] Skipping MMseqs2 for {path}: File empty or missing.")
    
    return valid_paths




localrules:
    all_mmseqs2,

rule all_mmseqs2:
    input:
        get_all_gene_mmseqs2_annotation,
    output:
        touch("genomes/annotations/{dataset}/genes/mmseqs2_easy_search/finished"),
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),


