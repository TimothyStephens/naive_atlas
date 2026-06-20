import os


#######################
####               ####
####    EGG NOG    ####
####               ####
#######################

rule gene_eggNOG_mapper:
    input:
        eggnog_db_files=rules.download_eggNOG_files.output.files,
        faa="genomes/genes/{dataset}/{genome}.faa",
    output:
        seed=temp(
            "genomes/annotations/{dataset}/genes/{genome}.emapper.seed_orthologs"
        ),
        hits=temp(
            "genomes/annotations/{dataset}/genes/{genome}.emapper.hits"
        ),
        annot=temp(
            "genomes/annotations/{dataset}/genes/{genome}.emapper.annotations"
        ),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else rules.download_eggNOG_files.output.dir
        ),
        prefix=lambda wc, output: output.annot.replace(".emapper.annotations", ""),
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: lambda wc: get_resource(wc, None, 1, "gene_annot_eggnog", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "gene_annot_eggnog", "account"),
    container:
        "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1"
    benchmark:
        "benchmarks/genomes/annotations/{dataset}/genes/{genome}.emapper_homology_search_diamond.tsv",
    log:
        "logs/genomes/annotations/{dataset}/genes/{genome}.emapper_homology_search_diamond.log",
    shell:
        """
        (
        if [ {params.copyto_shm} == "t" ] ; then
            # Check if the files exist before copying
            if [ ! -e "{params.data_dir}/eggnog.db" ]; then
                cp {EGGNOG_DIR}/eggnog.db {params.data_dir}/eggnog.db
            else
                echo "File {params.data_dir}/eggnog.db already exists. Skipping copy."
            fi
            
            if [ ! -e "{params.data_dir}/eggnog_proteins.dmnd" ]; then
                cp {EGGNOG_DIR}/eggnog_proteins.dmnd {params.data_dir}/eggnog_proteins.dmnd
            else
                echo "File {params.data_dir}/eggnog_proteins.dmnd already exists. Skipping copy."
            fi
        fi
        
        emapper.py \\
            -m diamond \\
            --no_file_comments \\
            --data_dir {params.data_dir} \\
            --dbmem \\
            --override \\
            -o {params.prefix} \\
            --cpu {threads} \\
            --data_dir {params.data_dir}
        ) 1>{log} 2>&1
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
            valid_paths.append(expand(rules.gene_eggNOG_mapper.output.annot, dataset=wildcards.dataset, genome=genome)[0])
        else:
            logger.info(f"[DEBUG] Skipping eggNOG for {path}: File empty or missing.")
    
    return valid_paths


rule combine_gene_egg_nog_annotations:
    input:
        get_all_gene_eggnog,
    output:
        parquet="genomes/annotations/{dataset}/genes/emapper.parquet",
        tsv="genomes/annotations/{dataset}/genes/emapper.tsv.gz",
    log:
        "logs/genomes/annotations/{dataset}/genes/emapper_combine.log",
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

            combined.columns = [
                "Query",
                "Seed",
                "Seed_evalue",
                "Seed_Score",
                "eggNOG",
                "max_annot_lvl",
                "COG_cat",
                "Description",
                "Name",
                "GO_terms",
                "EC",
                "KO",
                "KEGG_Pathway",
                "KEGG_Module",
                "KEGG_Reaction",
                "KEGG_rclass",
                "BRITE",
                "KEGG_TC",
                "CAZy",
                "BiGG_Reaction",
                "PFAMs",
            ]
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
        results=temp("genomes/annotations/{dataset}/genes/{genome}.faa.mmseqs2_easy_search_{database_name}.m4.gz"),
        tmp=temp(directory("genomes/annotations/{dataset}/genes/{genome}.faa.mmseqs2_easy_search_{database_name}.tmp")),
    params:
        opts=config["mmseqs2_easy_search_opts"],
        results="genomes/annotations/{dataset}/genes/{genome}.faa.mmseqs2_easy_search_{database_name}.m4",
        format_output=config["mmseqs2_easy_search_format_output"],
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
        "logs/genomes/annotations/{dataset}/genes/mmseqs2_easy_search_{database_name}/{genome}.log",
    benchmark:
        "benchmarks/genomes/annotations/{dataset}/genes/mmseqs2_easy_search_{database_name}/{genome}.tsv",
    shell:
        """
        (
        /usr/local/bin/entrypoint easy-search \\
            --threads {threads} \\
            --split-memory-limit {resources.mem_gb}G \\
            --format-mode 4 \\
            --format-output {params.format_output} \\
            {params.opts} \\
            {input.faa} \\
            {input.database} \\
            {params.results} \\
            {output.tmp} \\
          && gzip -9 {params.results}
        ) 1>{log} 2>&1
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
            logger.info(f"[DEBUG] Skipping MMseqs2 for {path}: File empty or missing.")
    
    return valid_paths




localrules:
    mmseqs2_combine,

rule mmseqs2_combine:
    input:
        get_all_gene_mmseqs2_annotation,
    output:
        "genomes/annotations/{dataset}/genes/mmseqs2_easy_search_{database_name}.m4.gz",
    log:
        "logs/genomes/annotations/{dataset}/genes/mmseqs2_combine_{database_name}.log",
    params:
        format_output=config["mmseqs2_easy_search_format_output"],
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
                pd.read_csv(file, index_col=None, header=None, sep="\t")
                for file in input
            ]

            combined = pd.concat(Tables, axis=0)

            del Tables

            combined.columns = params['format_output'].split(',')
            combined = combined.astype(str)

            combined.to_csv(output[0], sep='\t', index=False)
        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


