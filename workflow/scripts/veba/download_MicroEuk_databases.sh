#!/bin/bash
# __version__ = "2024.9.21"
# MICROEUKAYROTIC_DATABASE_VERSION = "MicroEuk_v3"
# usage: bash download_MicroEuk_databases.sh /path/to/veba_database_destination/

# Create database
DATABASE_DIRECTORY=${1:-"."}
REALPATH_DATABASE_DIRECTORY=$(realpath $DATABASE_DIRECTORY)
SCRIPT_DIRECTORY=$(dirname "$0")

MAXIMUM_NUMBER_OF_CPU=$(python -c "from multiprocessing import cpu_count; print(cpu_count())")
N_JOBS=${3:-${MAXIMUM_NUMBER_OF_CPU}}

# Database structure
echo ". .. ... ..... ........ ............."
echo "Creating directories"
echo ". .. ... ..... ........ ............."

mkdir -vp $DATABASE_DIRECTORY

# Versions
DATE=$(date)
echo $DATE > ${DATABASE_DIRECTORY}/ACCESS_DATE

echo ". .. ... ..... ........ ............."
echo " * Processing Microeukaryotic MMSEQS2 database"
echo ". .. ... ..... ........ ............."

# Download MicroEuk_v3 from Zenodo
wget -v -O ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz https://zenodo.org/records/10139451/files/MicroEuk_v3.tar.gz?download=1 
tar xvzf ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz -C ${DATABASE_DIRECTORY}

# Source Taxonomy
cp -rf ${DATABASE_DIRECTORY}/MicroEuk_v3/source_taxonomy.tsv.gz ${DATABASE_DIRECTORY}/

# MicroEuk100
gzip -d ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa.gz
mmseqs createdb --compressed 1 ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa ${DATABASE_DIRECTORY}/MicroEuk100

# MicroEuk100.eukaryota_odb10
gzip -d ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.eukaryota_odb10.list.gz
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.eukaryota_odb10.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa \
    | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/MicroEuk100.eukaryota_odb10

# MicroEuk90
gzip -d -c ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90_clusters.tsv.gz \
    | cut -f1 | sort -u > ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90.list
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk90.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa \
    | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/MicroEuk90

# MicroEuk50
gzip -d -c ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50_clusters.tsv.gz \
    | cut -f1 | sort -u > ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50.list
seqkit grep -f ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk50.list ${DATABASE_DIRECTORY}/MicroEuk_v3/MicroEuk100.faa \
    | mmseqs createdb --compressed 1 stdin ${DATABASE_DIRECTORY}/MicroEuk50

# source_to_lineage.dict.pkl.gz
build_source_to_lineage_dictionary.py -i ${DATABASE_DIRECTORY}/MicroEuk_v3/source_taxonomy.tsv.gz -o ${DATABASE_DIRECTORY}/source_to_lineage.dict.pkl.gz

# target_to_source.dict.pkl.gz
build_target_to_source_dictionary.py -i ${DATABASE_DIRECTORY}/MicroEuk_v3/identifier_mapping.proteins.tsv.gz -o ${DATABASE_DIRECTORY}/target_to_source.dict.pkl.gz

# Remove intermediate files
rm -rf ${DATABASE_DIRECTORY}/MicroEuk_v3
rm -rf ${DATABASE_DIRECTORY}/MicroEuk_v3.tar.gz

echo -e " _    _ _______ ______  _______\n  \  /  |______ |_____] |_____|\n   \/   |______ |_____] |     |"
echo -e "...................................................."
echo -e "             Database Download Complete             "
echo -e "...................................................."
