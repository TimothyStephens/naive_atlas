# Dockerfile to containerize naiveATLAS workflow
FROM condaforge/mambaforge:latest


## (1/5) Set environment variables
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


## (2/5) Install Singulairty
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      software-properties-common \
      gnupg \
      ca-certificates \
    && add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && apt-get install -y apptainer \
    && rm -rf /var/lib/apt/lists/*


## (3/5) Create and set workflow working directory
WORKDIR /app
COPY . .
RUN sed -e 's@timothystephens/naive_atlas-envs:latest@timothystephens/naive_atlas-envs:v0.3.4-dirty@' -i workflow/Snakefile
RUN mamba env create --prefix /conda-envs/naive_atlas --file naive_atlasenv_cluster.yml -v
SHELL ["conda", "run", "-p", "/conda-envs/naive_atlas", "/bin/bash", "-c"]
RUN /conda-envs/naive_atlas/bin/pip install --prefix /conda-envs/naive_atlas --editable .


## (4/5) Verify the installation by checking the version
RUN naive_atlas --help


## (5/5) Set the default command to show MaxBin2 help
CMD ["naive_atlas", "--help"]

