import hashlib
import os


# this values are incuded in the snakefile
DBDIR = os.path.realpath(config["database_dir"])
CHECKMDIR = os.path.join(DBDIR, "checkm")
CHECKM_ARCHIVE = "checkm_data_v1.0.9.tar.gz"
CAT_DIR = os.path.join(DBDIR, "CAT")
CAT_flag_downloaded = os.path.join(CAT_DIR, "downloaded")

ZENODO_ARCHIVE = "1134890"
EGGNOG_VERSION = "5"
EGGNOG_DIR = os.path.join(DBDIR, "EggNOG_V" + EGGNOG_VERSION)

GTDB_VERSION = "R220"
GTDB_DATA_URL = "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
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
    CHECKM_ARCHIVE: "631012fa598c43fdeb88c619ad282c4d",
}


CHECKMFILES = [
    "%s/taxon_marker_sets.tsv" % CHECKMDIR,
    "%s/selected_marker_sets.tsv" % CHECKMDIR,
    "%s/pfam/tigrfam2pfam.tsv" % CHECKMDIR,
    "%s/pfam/Pfam-A.hmm.dat" % CHECKMDIR,
    "%s/img/img_metadata.tsv" % CHECKMDIR,
    "%s/hmms_ssu/SSU_euk.hmm" % CHECKMDIR,
    "%s/hmms_ssu/SSU_bacteria.hmm" % CHECKMDIR,
    "%s/hmms_ssu/SSU_archaea.hmm" % CHECKMDIR,
    "%s/hmms_ssu/createHMMs.py" % CHECKMDIR,
    "%s/hmms/phylo.hmm.ssi" % CHECKMDIR,
    "%s/hmms/phylo.hmm" % CHECKMDIR,
    "%s/hmms/checkm.hmm.ssi" % CHECKMDIR,
    "%s/hmms/checkm.hmm" % CHECKMDIR,
    "%s/genome_tree/missing_duplicate_genes_97.tsv" % CHECKMDIR,
    "%s/genome_tree/missing_duplicate_genes_50.tsv" % CHECKMDIR,
    "%s/genome_tree/genome_tree.taxonomy.tsv" % CHECKMDIR,
    "%s/genome_tree/genome_tree_reduced.refpkg/phylo_modelJqWx6_.json" % CHECKMDIR,
    "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.tre" % CHECKMDIR,
    "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.log" % CHECKMDIR,
    "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.fasta" % CHECKMDIR,
    "%s/genome_tree/genome_tree_reduced.refpkg/CONTENTS.json" % CHECKMDIR,
    "%s/genome_tree/genome_tree.metadata.tsv" % CHECKMDIR,
    "%s/genome_tree/genome_tree_full.refpkg/phylo_modelEcOyPk.json" % CHECKMDIR,
    "%s/genome_tree/genome_tree_full.refpkg/genome_tree.tre" % CHECKMDIR,
    "%s/genome_tree/genome_tree_full.refpkg/genome_tree.log" % CHECKMDIR,
    "%s/genome_tree/genome_tree_full.refpkg/genome_tree.fasta" % CHECKMDIR,
    "%s/genome_tree/genome_tree_full.refpkg/CONTENTS.json" % CHECKMDIR,
    "%s/genome_tree/genome_tree.derep.txt" % CHECKMDIR,
    "%s/.dmanifest" % CHECKMDIR,
    "%s/distributions/td_dist.txt" % CHECKMDIR,
    "%s/distributions/gc_dist.txt" % CHECKMDIR,
    "%s/distributions/cd_dist.txt" % CHECKMDIR,
]


def get_eggnog_db_file():
    return ancient(
        expand(
            "{path}/{files}",
            path=EGGNOG_DIR,
            files=["eggnog.db", "eggnog_proteins.dmnd"],
        )
    )


localrules:
    download,
    download_eggNOG_files,
    download_atlas_files,
    download_checkm_data,


ruleorder: download_eggNOG_files > download_atlas_files


rule download:
    input:
        expand(
            "{dir}/{filename}", dir=DBDIR, filename=["adapters.fa", "phiX174_virus.fa"]
        ),
        get_eggnog_db_file(),
        f"{DBDIR}/CheckM2",
        os.path.join(GTDBTK_DATA_PATH, "downloaded_success"),


rule download_eggNOG_files:
    output:
        f"{EGGNOG_DIR}/eggnog.db",
        f"{EGGNOG_DIR}/eggnog_proteins.dmnd",
    conda:
        "../envs/eggNOG.yaml"
    shell:
        f"download_eggnog_data.py -yf --data_dir {EGGNOG_DIR} "


rule download_atlas_files:
    output:
        f"{DBDIR}/{{filename}}",
    wildcard_constraints:
        filename="[A-Za-z0-9_.]+",
    run:
        shell(
            "wget -O {output} 'https://zenodo.org/record/{ZENODO_ARCHIVE}/files/{wildcards.filename}' "
        )
        if not FILES[wildcards.filename] == md5(output[0]):
            raise OSError(2, "Invalid checksum", output[0])


rule download_checkm_data:
    output:
        tar=temp(CHECKM_ARCHIVE),
        files=CHECKMFILES,
    params:
        path=CHECKMDIR,
    run:
        shell(
            "wget -O {output.tar} 'https://zenodo.org/record/{ZENODO_ARCHIVE}/files/{CHECKM_ARCHIVE}' "
        )
        if not FILES[CHECKM_ARCHIVE] == md5(output.tar):
            raise OSError(2, "Invalid checksum", CHECKM_ARCHIVE)

        shell("tar -zxf {output.tar} --directory {params.path}")


localrules:
    initialize_checkm,


rule initialize_checkm:
    input:
        ancient(CHECKMFILES),
    output:
        touched_output=touch("logs/checkm_init.txt"),
    params:
        database_dir=CHECKMDIR,
    conda:
        "../envs/checkm.yaml"
    log:
        "logs/initialize_checkm.log",
    shell:
        "checkm data setRoot {params.database_dir} &> {log} "


localrules:
    download_gtdb,


rule download_gtdb:
    output:
        temp(f"{GTDBTK_DATA_PATH}/gtdb_data.tar.gz"),
    conda:
        "../envs/gtdbtk.yaml"
    resources:
        time=config["simplejob_runtime"],
    log:
        "logs/download/gtdbtk.log",
    shell:
        "wget --no-check-certificate {GTDB_DATA_URL} -O {output} &> {log} "


rule extract_gtdb:
    input:
        rules.download_gtdb.output,
    output:
        touch(os.path.join(GTDBTK_DATA_PATH, "downloaded_success")),
    conda:
        "../envs/gtdbtk.yaml"
    resources:
        time=config["simplejob_runtime"],
    log:
        "logs/download/gtdbtk_untar.log",
    shell:
        'tar -xzvf {input} -C "{GTDBTK_DATA_PATH}" --strip 1 2> {log}; '


rule checkm2_download_db:
    output:
        directory(f"{DBDIR}/CheckM2"),
    conda:
        "../envs/checkm2.yaml"
    log:
        "logs/download/checkm2.log",
    resources:
        time=config["simplejob_runtime"],
    shell:
        " checkm2 database --download --path {output} "
        " &>> {log}"


localrules:
    veba_download,

rule veba_download:
    output:
        dbdir=directory(f"{DBDIR}/veba"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    log:
        "logs/binning/download_veba_binning_database.log",
    benchmark:
        "logs/benchmarks/binning/download_veba_binning_database.tsv"
    conda:
        "../envs/VEBA-database_env.yml"
    shell:
        f"bash {workflow_folder}/scripts/veba/download_databases.sh"
        " {output.dbdir}"
        " {threads}"
        " &> {log}"


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
