import hashlib
import os


# This values are incuded in the snakefile
DBDIR = os.path.realpath(config["database_dir"])

GTDB_DATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
GTDBTK_DATA_PATH = os.path.join(DBDIR, "GTDB_R232")



localrules:
    all_downloads,

rule all_downloads:
    input:
        # Binning
        f"{DBDIR}/CheckM2",
        f"{DBDIR}/MDMcleaner",
        f"{DBDIR}/busco_lineages",
        f"{DBDIR}/geNomad",
        # Annotation
        f"{DBDIR}/EggNOG",
        os.path.join(f"{DBDIR}/MMseqs2", config["mmseqs2_database_name"]),
        os.path.join(GTDBTK_DATA_PATH, "downloaded_success"),
        # Gene Prediction
        f"{DBDIR}/MicroEuk",
        f"{DBDIR}/bakta/db",
    output:
        touch(f"{DBDIR}/finished")
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),




###############################
####                       ####
####        Binning        ####
####                       ####
###############################

rule checkm2_download_db:
    output:
        dbdir=directory(f"{DBDIR}/CheckM2"),
    conda:
        "../envs/checkm2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/checkm2.log",
    benchmark:
        "benchmarks/download/checkm2.tsv",
    shell:
        """
        (checkm2 database --download --path {output} --no_write_json_db) 1>{log} 2>&1
        """


rule mdmcleaner_download_db:
    output:
        dbdir=directory(f"{DBDIR}/MDMcleaner"),
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/mdmcleaner_database.log",
    benchmark:
        "benchmarks/download/mdmcleaner_database.tsv",
    container:
        "docker://timothystephens/mdmcleaner:0.8.7-TGSv3",
    shell:
        """
        (mdmcleaner makedb --outdir {output.dbdir}) 1>{log} 2>&1
        """


rule busco_download_db:
    output:
        dbdir=directory(f"{DBDIR}/busco_lineages"),
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/busco_lineages.log",
    benchmark:
        "benchmarks/download/busco_lineages.tsv",
    container:
        "docker://timothystephens/busco:6.0.0-TGSv1",
    shell:
        """
        (	
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export PYTHONPATH="$CONDA_PREFIX/lib/python3.7/site-packages"
        busco -q --download_path {output} --download all
        ) 1>{log} 2>&1
        """


rule genomad_download_db:
    output:
        dbdir=directory(f"{DBDIR}/geNomad"),
    params:
        db_version="v1.2",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/genomad_lineages.log",
    benchmark:
        "benchmarks/download/genomad_lineages.tsv",
    container:
        "docker://antoniopcamargo/genomad:1.11.0",
    shell:
        """
        (
        export PATH="/opt/conda/bin:$PATH"
        mkdir -p {output.dbdir}
        genomad download-database {output.dbdir}
        ) 1>{log} 2>&1
        """



###############################
####                       ####
####       Annotation      ####
####                       ####
###############################

rule download_eggNOG_files:
    output:
        files=[f"{DBDIR}/EggNOG/eggnog.db", f"{DBDIR}/EggNOG/eggnog_proteins.dmnd"],
        dir=f"{DBDIR}/EggNOG",
    params:
        eggnog_dir=f"{DBDIR}/EggNOG",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/download_eggNOG_files.log",
    benchmark:
        "benchmarks/download/download_eggNOG_files.tsv",
    container:
        "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1"
    shell:
        """
        (download_eggnog_data.py -yf --data_dir {params.eggnog_dir}) 1>{log} 2>&1
        """


rule gtdb_download_db:
    output:
        temp(f"{GTDBTK_DATA_PATH}/gtdb_data.tar.gz"),
    params:
        gtdb_data_url=f"{GTDB_DATA_URL}",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/gtdbtk.log",
    benchmark:
        "benchmarks/download/gtdbtk.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        (wget --no-check-certificate {params.gtdb_data_url} -O {output}) 1>{log} 2>&1
        """


rule gtdb_extract:
    input:
        rules.gtdb_download_db.output,
    output:
        touch(os.path.join(GTDBTK_DATA_PATH, "downloaded_success")),
    conda:
        "../envs/gtdbtk.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/gtdbtk_untar.log",
    benchmark:
        "benchmarks/download/gtdbtk_untar.tsv",
    shell:
        """
        (tar -xzvf {input} -C "{GTDBTK_DATA_PATH}" --strip 1) 1>{log} 2>&1
        """


rule mmseqs2_download:
    output:
        dbdir=directory(f"{DBDIR}/MMseqs2"),
        database=os.path.join(f"{DBDIR}/MMseqs2", config["mmseqs2_database_name"]),
    params:
        mmseqs2_database=config["mmseqs2_database"],
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/download_MMseqs2_database.log",
    benchmark:
        "benchmarks/download/download_MetaEuk_database.tsv",
    container:
        # Need a specific version of mmseqs2 other wise easy-taxonomy fails.
        "docker://timothystephens/mmseqs2:113e3212c137d026e297c7540e1fcd039f6812b1_rev1"
    shell:
        """
        (
        export TMPDIR='{output.dbdir}/tmp'; 
        mmseqs databases {params.mmseqs2_database} {output.database} {output.dbdir}/tmp \\
            --compressed 1 \\
            --threads {threads} \\
        && rm -fr {output.dbdir}/tmp
        ) 1>{log} 2>&1
        """



###############################
####                       ####
####    Gene Prediction    ####
####                       ####
###############################

rule microeukaryotic_mmseqs2_db:
    output:
        dbdir=directory(f"{DBDIR}/MicroEuk"),
    params:
        workflow_folder=os.path.dirname(os.path.abspath(workflow.snakefile)),
    log:
        "logs/download/microeukaryotic_mmseqs2.log",
    benchmark:
        "benchmarks/download/microeukaryotic_mmseqs2.tsv",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    conda:
        "../envs/MicroEuk.yaml"
    shell:
        """
        ({params.workflow_folder}/../scripts/veba/download_MicroEuk_databases.sh {output.dbdir}) 1>{log} 2>&1
        """


rule bakta_download_db:
    output:
        dbdir=directory(f"{DBDIR}/bakta/db"),
    params:
        wd=f"{DBDIR}/bakta",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/bakta.log",
    benchmark:
        "benchmarks/download/bakta.tsv",
    conda:
        "../envs/gene_prediction_bacteria.yaml"
    shell:
        """
        (
        mkdir -p {params.wd}; cd {params.wd}/
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        bakta_db download --type full
        ) 1>{log} 2>&1
        """



onsuccess:
    print("All databases have downloaded and validated successfully")


onerror:
    print("An error occurred while downloading reference databases.")
    print(
        "ATLAS databases can be manually downloaded from: https://zenodo.org/record/%s"
        % ZENODO_ARCHIVE
    )
    print(
        "eggNOG databases can be manually downloaded from: http://eggnogdb.embl.de/download/emapperdb-%s"
        % EGGNOG_VERSION
    )
    print(
        "CAT databases can be manually downloaded from: https://github.com/dutilh/CAT"
    )
