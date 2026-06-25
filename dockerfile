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


## (3/6) Create and set workflow working directory
WORKDIR /app
COPY . .
RUN sed -e 's@timothystephens/naive_atlas-envs:latest@timothystephens/naive_atlas-envs:v0.3.4-35-g03cb9c8-dirty@' -i workflow/Snakefile
RUN mamba env create --prefix /conda-envs/naive_atlas --file naive_atlasenv_cluster.yml -v
SHELL ["conda", "run", "-p", "/conda-envs/naive_atlas", "/bin/bash", "-c"]
RUN /conda-envs/naive_atlas/bin/pip install --prefix /conda-envs/naive_atlas --editable .


## (4/6) Verify the installation by checking the version
RUN naive_atlas --help


## (5/6) Set the default command to test naive_atlas
CMD ["naive_atlas", "--help"]

## (6/6) Download containers
RUN mkdir -p /singulairty \
    && apptainer pull --dir /singulairty 8495bd0424646c72faf1123bc3b89544.sif "docker://antoniopcamargo/genomad:1.11.0" \
    && apptainer pull --dir /singulairty 622de81f72bce099ff2a1b9d977b3801.sif "docker://chrishah/cdhit:v4.8.1" \
    && apptainer pull --dir /singulairty a537d2092becf2c9685c369e295c91f4.sif "docker://ghcr.io/soedinglab/metaeuk:7-bba0d80" \
    && apptainer pull --dir /singulairty 6008b3d8e1ae9bb57e31b4426fe90284.sif "docker://ghcr.io/soedinglab/mmseqs2:18-8cc5c" \
    && apptainer pull --dir /singulairty 35740ed4f0505cb05d27492add3564b1.sif "docker://mambaorg/micromamba:2.5.0-cuda11.8.0-ubuntu22.04" \
    && apptainer pull --dir /singulairty c3d1e9552dee486c6b84387682fe5ad3.sif "docker://timothystephens/busco:6.1.0-TGSv1" \
    && apptainer pull --dir /singulairty d75eff3ead4a99ccfabdb60d3bbc8169.sif "docker://timothystephens/eggnog-mapper:2.1.13-TGSv1" \
    && apptainer pull --dir /singulairty ad6661eef59642c09d33e3dda8781db4.sif "docker://timothystephens/maxbin2:2.2.7-TGSv5" \
    && apptainer pull --dir /singulairty c1b1e97d1c1b7d133d5779600411ef52.sif "docker://timothystephens/mdmcleaner:0.8.7-TGSv4" \
    && apptainer pull --dir /singulairty aad2074efde11488d9c164db909d1372.sif "docker://timothystephens/mmseqs2:113e3212c137d026e297c7540e1fcd039f6812b1_rev1" \
    && apptainer cache clean --force
RUN apptainer pull --dir /singulairty b0cbfcaa12aed20fbe0a5297d12b5d21.sif "docker-archive:///app/naive_atlas-envs.tar" \
    && rm -f /app/naive_atlas-envs.tar \
    && apptainer cache clean --force
