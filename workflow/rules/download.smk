import hashlib
import os


# this values are incuded in the snakefile
DBDIR = os.path.realpath(config["database_dir"])

ZENODO_ARCHIVE = "1134890"
EGGNOG_VERSION = "5"
EGGNOG_DIR = os.path.join(DBDIR, "EggNOG_V" + EGGNOG_VERSION)

GTDB_VERSION = "R226"
GTDB_DATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz"
#GTDB_DATA_URL = "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
GTDBTK_DATA_PATH = os.path.join(DBDIR, "GTDB_" + GTDB_VERSION)


def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# note: saving OG_fasta.tar.gz in order to not create secondary "success" file
FILES = {
    "adapters.fa": "ae839dc79cfb855a1b750a0d593fe01e",
    "phiX174_virus.fa": "82516880142e8c89b466bc6118696c47",
    "refseq.db": "42b8976656f2cfd661b8a299d6e24c19",
    "refseq.dmnd": "c01facc7e397270ccb796ea799a09108",
    "refseq.tree": "469fcbeb15dd0d4bf8f1677682bde157",
    "silva_rfam_all_rRNAs.fa": "f102e35d9f48eabeb0efe9058559bc66",
    "eggnog.db": "7923d3bb7eca8e0e8f122be4b5ca6997",
    "eggnog_proteins.dmnd": "64fefa838833a6f3e220a06fb9d403cd",
}


def get_eggnog_db_file():
    return ancient(
        expand(
            "{path}/{files}",
            path=EGGNOG_DIR,
            files=["eggnog.db", "eggnog_proteins.dmnd"],
        )
    )


ruleorder: download_eggNOG_files > download_atlas_files


localrules:
    all_downloads,

rule all_downloads:
    input:
        expand("{dir}/{filename}", dir=DBDIR, filename=["adapters.fa", "phiX174_virus.fa"]),
        # Binning
        f"{DBDIR}/CheckM2",
        f"{DBDIR}/MDMcleaner",
        f"{DBDIR}/busco_lineages",
        f"{DBDIR}/geNomad",
        # Annotation
        get_eggnog_db_file(),
        os.path.join(f"{DBDIR}/MMseqs2", config["mmseqs2_database_name"]),
        f"{DBDIR}/DRAM/db/",
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


rule download_atlas_files:
    output:
        f"{DBDIR}/{{filename}}",
    wildcard_constraints:
        filename="[A-Za-z0-9_.]+",
    log:
        "logs/download/download_atlas_file_{filename}.log",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    benchmark:
        "benchmarks/download/download_atlas_file_{filename}.tsv"
    run:
        shell(
            "wget -O {output} 'https://zenodo.org/record/{ZENODO_ARCHIVE}/files/{wildcards.filename}' "
        )
        if not FILES[wildcards.filename] == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])



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
        "benchmarks/download/checkm2.tsv"
    shell:
        """
        checkm2 database --download --path {output} --no_write_json_db &> {log}
        """


rule mdmcleaner_download_db:
    output:
        dbdir=directory(f"{DBDIR}/MDMcleaner"),
    log:
        "logs/download/mdmcleaner_database.log",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    benchmark:
        "benchmarks/download/mdmcleaner_database.tsv"
    container:
        "docker://timothystephens/mdmcleaner:0.8.7-TGSv2",
    #conda:
    #    "../envs/mdmcleaner.yaml"
    shell:
        """
        mdmcleaner makedb --outdir {output.dbdir} &> {log}
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
        "benchmarks/download/busco_lineages.tsv"
    conda:
        "../envs/busco.yaml"
    shell:
        """
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export PYTHONPATH="$CONDA_PREFIX/lib/python3.7/site-packages"
        busco -q --download_path {output} --download all &> {log}
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
        "benchmarks/download/genomad_lineages.tsv"
    #conda:
    #    "../envs/genomad.yaml"
    container:
        "docker://antoniopcamargo/genomad:1.11.0",
    shell:
        """
        (
        export PATH="/opt/conda/bin:$PATH"
        mkdir -p {output.dbdir}
        genomad download-database {output.dbdir}
        ) &> {log}
        """



###############################
####                       ####
####       Annotation      ####
####                       ####
###############################

rule download_eggNOG_files:
    output:
        f"{EGGNOG_DIR}/eggnog.db",
        f"{EGGNOG_DIR}/eggnog_proteins.dmnd",
    params:
        eggnog_dir=f"{EGGNOG_DIR}",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/download/download_eggNOG_files.log",
    benchmark:
        "benchmarks/download/download_eggNOG_files.tsv"
    #conda:
    #    "../envs/eggNOG.yaml"
    container:
        "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1"
    shell:
        """
        download_eggnog_data.py -yf --data_dir {params.eggnog_dir} &> {log}
        """


rule dram_download:
    output:
        dbdir=directory(f"{DBDIR}/DRAM/db/"),
        config=f"{DBDIR}/DRAM/DRAM.config",
    threads: lambda wc: get_resource(wc, None, 1, "download", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "download", "account"),
    log:
        "logs/dram/download_dram.log",
    benchmark:
        "benchmarks/dram/download_dram.tsv"
    conda:
        "../envs/dram.yaml"
    shell:
        " DRAM-setup.py prepare_databases "
        " --output_dir {output.dbdir} "
        " --threads {threads} "
        " --verbose "
        " --skip_uniref "
        " &> {log} "
        " ; "
        " DRAM-setup.py export_config --output_file {output.config}"


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
        "benchmarks/download/gtdbtk.tsv"
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        wget --no-check-certificate {params.gtdb_data_url} -O {output} &> {log}
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
        "benchmarks/download/gtdbtk_untar.tsv"
    shell:
        """
        tar -xzvf {input} -C "{GTDBTK_DATA_PATH}" --strip 1 &> {log}
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
        "benchmarks/download/download_MetaEuk_database.tsv"
    container:
        # Need a specific version of mmseqs2 other wise easy-taxonomy fails.
        "docker://timothystephens/mmseqs2:113e3212c137d026e297c7540e1fcd039f6812b1_rev1"
    shell:
        "mmseqs databases {params.mmseqs2_database} {output.database} {output.dbdir}/tmp "
        " --compressed 1 "
        " --threads {threads} "
        " &> {log}"



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
        "benchmarks/download/microeukaryotic_mmseqs2.tsv"
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
        {params.workflow_folder}/../scripts/veba/download_MicroEuk_databases.sh {output.dbdir} &> {log}
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
        "benchmarks/download/bakta.tsv"
    conda:
        "../envs/gene_prediction_bacteria.yaml"
    shell:
        """
        (
        mkdir -p {params.wd}; cd {params.wd}/
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        bakta_db download --type full
        ) &> {log}
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
