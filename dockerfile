# Dockerfile to containerize naiveATLAS workflow
FROM condaforge/mambaforge:latest


## (1/6) Set environment variables
# So we dont have to interactivly configure tzdata
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
ENV DEBIAN_FRONTEND=noninteractive

ENV CONDA_DIR=/conda-envs
ENV PATH=/conda-envs/naive_atlas/bin:$PATH
# Enhance mamba network robustness
ENV MAMBA_NO_LOW_SPEED_LIMIT=1
ENV CONDA_REMOTE_READ_TIMEOUT_SECS=300

RUN apt-get update && \
    apt-get install -y \
      zlib1g zlib1g-dev build-essential gcc && \
    rm -rf /var/lib/apt/lists/*


## (2/6) Install Singulairty
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      software-properties-common \
      gnupg \
      ca-certificates \
    && add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && apt-get install -y apptainer \
    && rm -rf /var/lib/apt/lists/*


## (3/6) Install each workflow package
# Conda environment:
#   source: workflow/envs/assembly.yaml
#   prefix: /conda-envs/2c8123cf913842181227d9a95d916f6c
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - seqkit=2.4.0
#       - spades=4.2.0
#       - megahit=1.2.9
#       - flye=2.9.6
#       - metamdbg=1.3.1
COPY workflow/envs/assembly.yaml /conda-envs/2c8123cf913842181227d9a95d916f6c.yaml
RUN mamba env create  \
      --prefix /conda-envs/2c8123cf913842181227d9a95d916f6c \
      --file   /conda-envs/2c8123cf913842181227d9a95d916f6c.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/busco.yaml
#   prefix: /conda-envs/43a1a2a96f2107f381edf1301001e5fa
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - busco=6.0.0
#   - biopython=1.79
COPY workflow/envs/busco.yaml /conda-envs/43a1a2a96f2107f381edf1301001e5fa.yaml
RUN mamba env create  \
      --prefix /conda-envs/43a1a2a96f2107f381edf1301001e5fa \
      --file   /conda-envs/43a1a2a96f2107f381edf1301001e5fa.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/checkm2.yaml
#   prefix: /conda-envs/df153c97c6b9f88f8d483af0a35e1402
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkm2>=1.1.0
#     - biopython=1.79
COPY workflow/envs/checkm2.yaml /conda-envs/df153c97c6b9f88f8d483af0a35e1402.yaml
RUN mamba env create  \
      --prefix /conda-envs/df153c97c6b9f88f8d483af0a35e1402 \
      --file   /conda-envs/df153c97c6b9f88f8d483af0a35e1402.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/checkm.yaml
#   prefix: /conda-envs/70c51d37c49a9da0979ed4257f1d558c
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - checkm-genome=1.2.4
COPY workflow/envs/checkm.yaml /conda-envs/70c51d37c49a9da0979ed4257f1d558c.yaml
RUN mamba env create  \
      --prefix /conda-envs/70c51d37c49a9da0979ed4257f1d558c \
      --file   /conda-envs/70c51d37c49a9da0979ed4257f1d558c.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/checkv.yaml
#   prefix: /conda-envs/daba6532d8e851dd42bafd749dd61eed
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - checkv=1.0.3
COPY workflow/envs/checkv.yaml /conda-envs/daba6532d8e851dd42bafd749dd61eed.yaml
RUN mamba env create  \
      --prefix /conda-envs/daba6532d8e851dd42bafd749dd61eed \
      --file   /conda-envs/daba6532d8e851dd42bafd749dd61eed.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/compleasm.yaml
#   prefix: /conda-envs/d1e883fbc33be95ce7f00cbb24af432b
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - compleasm=0.2.7
COPY workflow/envs/compleasm.yaml /conda-envs/d1e883fbc33be95ce7f00cbb24af432b.yaml
RUN mamba env create  \
      --prefix /conda-envs/d1e883fbc33be95ce7f00cbb24af432b \
      --file   /conda-envs/d1e883fbc33be95ce7f00cbb24af432b.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/coverm.yaml
#   prefix: /conda-envs/0608d7638ce1bc0c43d56ab74f3bfdf7
#   name: coverm
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - coverm=0.7.0
#     - samtools=1.16.1
COPY workflow/envs/coverm.yaml /conda-envs/0608d7638ce1bc0c43d56ab74f3bfdf7.yaml
RUN mamba env create  \
      --prefix /conda-envs/0608d7638ce1bc0c43d56ab74f3bfdf7 \
      --file   /conda-envs/0608d7638ce1bc0c43d56ab74f3bfdf7.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/dastool.yaml
#   prefix: /conda-envs/cdaccecb0d819c466e118c369559b513
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - das_tool=1.1.7
#       - python=3.11.4
COPY workflow/envs/dastool.yaml /conda-envs/cdaccecb0d819c466e118c369559b513.yaml
RUN mamba env create  \
      --prefix /conda-envs/cdaccecb0d819c466e118c369559b513 \
      --file   /conda-envs/cdaccecb0d819c466e118c369559b513.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/dram.yaml
#   prefix: /conda-envs/ca7d65e3e4c5c9a6e6b02e648cc3601b
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python >=3.8
#     - altair >=4
#     - networkx
#     - numpy
#     - openpyxl
#     - pandas >=1.5, <2
#     - scikit-bio >=0.5.8, <0.6
#     - sqlalchemy
#     - prodigal
#     - scipy >=1.9
#     - mmseqs2 >10.6d92c
#     - hmmer
#     - trnascan-se >=2
#     - barrnap
#     - ruby
#     - parallel
#     - wget
#     - curl
#     - pip
#     - pip:
#       - git+https://github.com/SilasK/DRAM.git
COPY workflow/envs/dram.yaml /conda-envs/ca7d65e3e4c5c9a6e6b02e648cc3601b.yaml
RUN mamba env create  \
      --prefix /conda-envs/ca7d65e3e4c5c9a6e6b02e648cc3601b \
      --file   /conda-envs/ca7d65e3e4c5c9a6e6b02e648cc3601b.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/eggNOG.yaml
#   prefix: /conda-envs/1b58c0a8c1fdb8a1c216a818a24d8134
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - eggnog-mapper=2.1.13
#       - python=3.11
#       - diamond =2.1
#       - wget # to download_eggnog_data on macOS
#       - pandas>=1.5,<2
COPY workflow/envs/eggNOG.yaml /conda-envs/1b58c0a8c1fdb8a1c216a818a24d8134.yaml
RUN mamba env create  \
      --prefix /conda-envs/1b58c0a8c1fdb8a1c216a818a24d8134 \
      --file   /conda-envs/1b58c0a8c1fdb8a1c216a818a24d8134.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/fasta.yaml
#   prefix: /conda-envs/95c31fa251d5138322098a3958f194bb
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - pyfastx=0.9
#       - pandas=1.2
#       - pyarrow
#       - biopython
COPY workflow/envs/fasta.yaml /conda-envs/95c31fa251d5138322098a3958f194bb.yaml
RUN mamba env create  \
      --prefix /conda-envs/95c31fa251d5138322098a3958f194bb \
      --file   /conda-envs/95c31fa251d5138322098a3958f194bb.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gene_prediction_archaea.yaml
#   prefix: /conda-envs/e3fe7432bc0701d381d18197703d4940
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#     - prokka=1.15.6
COPY workflow/envs/gene_prediction_archaea.yaml /conda-envs/e3fe7432bc0701d381d18197703d4940.yaml
RUN mamba env create  \
      --prefix /conda-envs/e3fe7432bc0701d381d18197703d4940 \
      --file   /conda-envs/e3fe7432bc0701d381d18197703d4940.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gene_prediction_bacteria.yaml
#   prefix: /conda-envs/9fe58690fca93f66575486c3d2fb6b65
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#     - bakta=1.12.0
COPY workflow/envs/gene_prediction_bacteria.yaml /conda-envs/9fe58690fca93f66575486c3d2fb6b65.yaml
RUN mamba env create  \
      --prefix /conda-envs/9fe58690fca93f66575486c3d2fb6b65 \
      --file   /conda-envs/9fe58690fca93f66575486c3d2fb6b65.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gene_prediction_eukaryotic.yaml
#   prefix: /conda-envs/6ea354a57160f7e4ab9be7d9e5e938e1
#   channels:
#       - jolespin
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#     - barrnap=0.9
#     - metaeuk=5.34c21f2
#     - pyrodigal=2.1.0
#     - seqkit=2.4.0
#     - tiara=1.0.3
#     - trnascan-se=2.0.9
#     - pandas
#     - genopype
COPY workflow/envs/gene_prediction_eukaryotic.yaml /conda-envs/6ea354a57160f7e4ab9be7d9e5e938e1.yaml
RUN mamba env create  \
      --prefix /conda-envs/6ea354a57160f7e4ab9be7d9e5e938e1 \
      --file   /conda-envs/6ea354a57160f7e4ab9be7d9e5e938e1.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gene_prediction_plasmid.yaml
#   prefix: /conda-envs/9fe58690fca93f66575486c3d2fb6b65
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#     - bakta=1.12.0
# NOTRUN: Seen md5sum previously

# Conda environment:
#   source: workflow/envs/gene_prediction_virus.yaml
#   prefix: /conda-envs/e3fe7432bc0701d381d18197703d4940
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#     - prokka=1.15.6
# NOTRUN: Seen md5sum previously

# Conda environment:
#   source: workflow/envs/genomad.yaml
#   prefix: /conda-envs/9999a65bd6e217a51e897724fe66d532
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - genomad=1.11.2
COPY workflow/envs/genomad.yaml /conda-envs/9999a65bd6e217a51e897724fe66d532.yaml
RUN mamba env create  \
      --prefix /conda-envs/9999a65bd6e217a51e897724fe66d532 \
      --file   /conda-envs/9999a65bd6e217a51e897724fe66d532.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gtdbtk.yaml
#   prefix: /conda-envs/697aed27f4f24eed9d0d588b8ebca7d1
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - gtdbtk=2.6.1
COPY workflow/envs/gtdbtk.yaml /conda-envs/697aed27f4f24eed9d0d588b8ebca7d1.yaml
RUN mamba env create  \
      --prefix /conda-envs/697aed27f4f24eed9d0d588b8ebca7d1 \
      --file   /conda-envs/697aed27f4f24eed9d0d588b8ebca7d1.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/gunc.yaml
#   prefix: /conda-envs/fb001d455000d4529f88fe6f1e60f557
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - gunc=1.0.6
#   - pandas=1.5.1
COPY workflow/envs/gunc.yaml /conda-envs/fb001d455000d4529f88fe6f1e60f557.yaml
RUN mamba env create  \
      --prefix /conda-envs/fb001d455000d4529f88fe6f1e60f557 \
      --file   /conda-envs/fb001d455000d4529f88fe6f1e60f557.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/instrain.yaml
#   prefix: /conda-envs/702bbf3c6194bc68e4ab49e93128b050
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - instrain =1.9.0
#   - samtools
COPY workflow/envs/instrain.yaml /conda-envs/702bbf3c6194bc68e4ab49e93128b050.yaml
RUN mamba env create  \
      --prefix /conda-envs/702bbf3c6194bc68e4ab49e93128b050 \
      --file   /conda-envs/702bbf3c6194bc68e4ab49e93128b050.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/maxbin2.yaml
#   prefix: /conda-envs/bf8a3027fd33b3db7ab33f1778d137d2
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - maxbin2=2.2.7
#       - python=3.11.4
#       - pandas=2.0.2
#       - tqdm=4.65.0
COPY workflow/envs/maxbin2.yaml /conda-envs/bf8a3027fd33b3db7ab33f1778d137d2.yaml
RUN mamba env create  \
      --prefix /conda-envs/bf8a3027fd33b3db7ab33f1778d137d2 \
      --file   /conda-envs/bf8a3027fd33b3db7ab33f1778d137d2.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/mdmcleaner.yaml
#   prefix: /conda-envs/6bb6509d829d3f2521ae4ee640d4d39c
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - mdmcleaner=0.8.7
COPY workflow/envs/mdmcleaner.yaml /conda-envs/6bb6509d829d3f2521ae4ee640d4d39c.yaml
RUN mamba env create  \
      --prefix /conda-envs/6bb6509d829d3f2521ae4ee640d4d39c \
      --file   /conda-envs/6bb6509d829d3f2521ae4ee640d4d39c.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/metabat2.yaml
#   prefix: /conda-envs/c8cd60a24df25faa448d5d1ef3ee71b4
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - metabat2 =2.15
#       - python=3.11.4
#       - pandas=2.0.2
#       - tqdm=4.65.0
#       - seqkit=2.4.0
COPY workflow/envs/metabat2.yaml /conda-envs/c8cd60a24df25faa448d5d1ef3ee71b4.yaml
RUN mamba env create  \
      --prefix /conda-envs/c8cd60a24df25faa448d5d1ef3ee71b4 \
      --file   /conda-envs/c8cd60a24df25faa448d5d1ef3ee71b4.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/metaeuk.yaml
#   prefix: /conda-envs/cf25e640d2e3a420390ff90aef97788b
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - metaeuk=7.bba0d80
COPY workflow/envs/metaeuk.yaml /conda-envs/cf25e640d2e3a420390ff90aef97788b.yaml
RUN mamba env create  \
      --prefix /conda-envs/cf25e640d2e3a420390ff90aef97788b \
      --file   /conda-envs/cf25e640d2e3a420390ff90aef97788b.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/MicroEuk.yaml
#   prefix: /conda-envs/ee0886e2171ef9ebda6c37117672b7ee
#   channels:
#     - conda-forge
#     - bioconda
#     - default
#   dependencies:
#     - mmseqs2=14.7e284
#     - seqkit=2.4.0
#     - tqdm
#     - pandas
COPY workflow/envs/MicroEuk.yaml /conda-envs/ee0886e2171ef9ebda6c37117672b7ee.yaml
RUN mamba env create  \
      --prefix /conda-envs/ee0886e2171ef9ebda6c37117672b7ee \
      --file   /conda-envs/ee0886e2171ef9ebda6c37117672b7ee.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/minimap.yaml
#   prefix: /conda-envs/d91311468e70edfc304cd7bf1b2da061
#   name: minimap2
#   channels:
#     - conda-forge
#     - bioconda
#     - default
#   dependencies:
#     - minimap2=2.30
#     - samtools=1.22.1
COPY workflow/envs/minimap.yaml /conda-envs/d91311468e70edfc304cd7bf1b2da061.yaml
RUN mamba env create  \
      --prefix /conda-envs/d91311468e70edfc304cd7bf1b2da061 \
      --file   /conda-envs/d91311468e70edfc304cd7bf1b2da061.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/mmseqs2.yaml
#   prefix: /conda-envs/6e92808fa97389422e191eaf3314f0d0
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - mmseqs2=17.b804f
#       - pigz
COPY workflow/envs/mmseqs2.yaml /conda-envs/6e92808fa97389422e191eaf3314f0d0.yaml
RUN mamba env create  \
      --prefix /conda-envs/6e92808fa97389422e191eaf3314f0d0 \
      --file   /conda-envs/6e92808fa97389422e191eaf3314f0d0.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/prodigal.yaml
#   prefix: /conda-envs/42d7cd011be38e39ac59b7a66e7e0ac4
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - prodigal =2.6.*
COPY workflow/envs/prodigal.yaml /conda-envs/42d7cd011be38e39ac59b7a66e7e0ac4.yaml
RUN mamba env create  \
      --prefix /conda-envs/42d7cd011be38e39ac59b7a66e7e0ac4 \
      --file   /conda-envs/42d7cd011be38e39ac59b7a66e7e0ac4.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/python.yaml
#   prefix: /conda-envs/317f09ac24cecfca6bde4b29eed8378b
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - python=3.12
#       - biopython
#       - pandas
#       - tqdm
COPY workflow/envs/python.yaml /conda-envs/317f09ac24cecfca6bde4b29eed8378b.yaml
RUN mamba env create  \
      --prefix /conda-envs/317f09ac24cecfca6bde4b29eed8378b \
      --file   /conda-envs/317f09ac24cecfca6bde4b29eed8378b.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/report.yaml
#   prefix: /conda-envs/6a803d63f4380786f8bba080530b7136
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - pandas=1.2
#       - plotly=5.3
COPY workflow/envs/report.yaml /conda-envs/6a803d63f4380786f8bba080530b7136.yaml
RUN mamba env create  \
      --prefix /conda-envs/6a803d63f4380786f8bba080530b7136 \
      --file   /conda-envs/6a803d63f4380786f8bba080530b7136.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/required_packages.yaml
#   prefix: /conda-envs/09db3890a6876fe231cc03c0fab2665c
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - python >=3.8, <3.11
#       - bbmap=38.87
#       - pigz
#       - bzip2
#       - pandas >=1.2, <2
#       - samtools
#       - sambamba
COPY workflow/envs/required_packages.yaml /conda-envs/09db3890a6876fe231cc03c0fab2665c.yaml
RUN mamba env create  \
      --prefix /conda-envs/09db3890a6876fe231cc03c0fab2665c \
      --file   /conda-envs/09db3890a6876fe231cc03c0fab2665c.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/seqkit.yaml
#   prefix: /conda-envs/d8dc31c42d25b37e80d9985723c6ae00
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - python=3.12
#       - seqkit=2.4.0
#       - pandas
COPY workflow/envs/seqkit.yaml /conda-envs/d8dc31c42d25b37e80d9985723c6ae00.yaml
RUN mamba env create  \
      --prefix /conda-envs/d8dc31c42d25b37e80d9985723c6ae00 \
      --file   /conda-envs/d8dc31c42d25b37e80d9985723c6ae00.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/skani.yaml
#   prefix: /conda-envs/126486bcd9fd1ab264dd80a64fb0db83
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - skani=0.3.1
COPY workflow/envs/skani.yaml /conda-envs/126486bcd9fd1ab264dd80a64fb0db83.yaml
RUN mamba env create  \
      --prefix /conda-envs/126486bcd9fd1ab264dd80a64fb0db83 \
      --file   /conda-envs/126486bcd9fd1ab264dd80a64fb0db83.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/species_clustering.yaml
#   prefix: /conda-envs/b2d08303db669404a1248cce0d48f843
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.11
#     - pandas=2
#     - pyarrow=11
#     - networkx=3.1
#     - scipy=1.10
COPY workflow/envs/species_clustering.yaml /conda-envs/b2d08303db669404a1248cce0d48f843.yaml
RUN mamba env create  \
      --prefix /conda-envs/b2d08303db669404a1248cce0d48f843 \
      --file   /conda-envs/b2d08303db669404a1248cce0d48f843.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/sra.yaml
#   prefix: /conda-envs/5f9497eb27bd1404656d539a10c472d6
#   channels:
#   - defaults
#   - bioconda
#   - conda-forge
#   dependencies:
#   - sra-tools
#   - pigz
#   - parallel-fastq-dump
COPY workflow/envs/sra.yaml /conda-envs/5f9497eb27bd1404656d539a10c472d6.yaml
RUN mamba env create  \
      --prefix /conda-envs/5f9497eb27bd1404656d539a10c472d6 \
      --file   /conda-envs/5f9497eb27bd1404656d539a10c472d6.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/tree.yaml
#   prefix: /conda-envs/8e639eadbcd374d0d9db53e61659d15b
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - python=3.12
#       - ete3=3.1.2
COPY workflow/envs/tree.yaml /conda-envs/8e639eadbcd374d0d9db53e61659d15b.yaml
RUN mamba env create  \
      --prefix /conda-envs/8e639eadbcd374d0d9db53e61659d15b \
      --file   /conda-envs/8e639eadbcd374d0d9db53e61659d15b.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/virus_filter.yaml
#   prefix: /conda-envs/f4fc1e04c6421bfd568c068aa4bf42b7
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - python=3.12
#       - seqkit=2.4.0
#       - biopython
#       - pandas
#       - tqdm
COPY workflow/envs/virus_filter.yaml /conda-envs/f4fc1e04c6421bfd568c068aa4bf42b7.yaml
RUN mamba env create  \
      --prefix /conda-envs/f4fc1e04c6421bfd568c068aa4bf42b7 \
      --file   /conda-envs/f4fc1e04c6421bfd568c068aa4bf42b7.yaml -vvv && \
    mamba clean -afy

# Conda environment:
#   source: workflow/envs/whokaryote.yaml
#   prefix: /conda-envs/2d2e6a31d6cda87ec1e7252fdbbda64e
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - whokaryote=1.1.2
COPY workflow/envs/whokaryote.yaml /conda-envs/2d2e6a31d6cda87ec1e7252fdbbda64e.yaml
RUN mamba env create  \
      --prefix /conda-envs/2d2e6a31d6cda87ec1e7252fdbbda64e \
      --file   /conda-envs/2d2e6a31d6cda87ec1e7252fdbbda64e.yaml -vvv && \
    mamba clean -afy


## (4/6) Create and set workflow working directory
WORKDIR /app
COPY . .
RUN mamba env create --prefix /conda-envs/naive_atlas --file naive_atlasenv_cluster.yml -v
SHELL ["conda", "run", "-p", "/conda-envs/naive_atlas", "/bin/bash", "-c"]
RUN /conda-envs/naive_atlas/bin/pip install --prefix /conda-envs/naive_atlas --editable .


## (5/6) Verify the installation by checking the version
RUN naive_atlas --help


## (6/6) Set the default command to show MaxBin2 help
CMD ["naive_atlas", "--help"]

