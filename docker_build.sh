####
#### Script to automatically containerize and build naiveATLAS into docker
####
set -euo pipefail


## Env vars
USER="timothystephens"
DOCKER=$(which docker)
DOCKERFILE="dockerfile"
VERSION=$(git describe --tags --dirty --always)
DEBUG="-v"

## Build dockerfile
echo -e "# Dockerfile to containerize naiveATLAS workflow
FROM condaforge/mambaforge:latest


## (1/6) Set environment variables
# So we dont have to interactivly configure tzdata
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
ENV DEBIAN_FRONTEND=noninteractive

ENV CONDA_DIR=/conda-envs
ENV PATH=/conda-envs/naive_atlas/bin:\$PATH
# Enhance mamba network robustness
ENV MAMBA_NO_LOW_SPEED_LIMIT=1
ENV CONDA_REMOTE_READ_TIMEOUT_SECS=300

RUN apt-get update && \\
    apt-get install -y \\
      zlib1g zlib1g-dev build-essential gcc && \\
    rm -rf /var/lib/apt/lists/*


## (2/6) Install Singulairty
RUN apt-get update && \\
    apt-get install -y --no-install-recommends \\
      software-properties-common \\
      gnupg \\
      ca-certificates \\
    && add-apt-repository -y ppa:apptainer/ppa \\
    && apt-get update \\
    && apt-get install -y apptainer \\
    && rm -rf /var/lib/apt/lists/*


## (3/6) Install each workflow package" > "$DOCKERFILE"

RUN=""
tmp_seen=$(mktemp)
for YAML in workflow/envs/*.yaml;
do

MD5SUM=$(md5sum "$YAML" | awk '{print $1}')

# Check if we have seen this md5sum previously
if ! grep -qFx "$MD5SUM" "$tmp_seen"; then
  # Not seen
  echo "# Conda environment:" >> "$DOCKERFILE"
  echo "#   source: $YAML" >> "$DOCKERFILE"
  echo "#   prefix: /conda-envs/$MD5SUM" >> "$DOCKERFILE"
  awk '{print "#   "$0}' "$YAML"   >> "$DOCKERFILE"
  echo -e "COPY $YAML /conda-envs/$MD5SUM.yaml" >> "$DOCKERFILE"
  echo -e "RUN mamba env create ${DEBUG} \\" >> "$DOCKERFILE"
  echo -e "      --prefix /conda-envs/$MD5SUM \\" >> "$DOCKERFILE"
  echo -e "      --file   /conda-envs/$MD5SUM.yaml -vvv && \\" >> "$DOCKERFILE"
  echo -e "    mamba clean -afy" >> "$DOCKERFILE"
  echo "" >> "$DOCKERFILE"
  
  echo "$MD5SUM" >> "$tmp_seen"
  
else
  # Seen
  echo "# Conda environment:" >> "$DOCKERFILE"
  echo "#   source: $YAML" >> "$DOCKERFILE"
  echo "#   prefix: /conda-envs/$MD5SUM" >> "$DOCKERFILE"
  awk '{print "#   "$0}' "$YAML"   >> "$DOCKERFILE"
  echo -e "# NOTRUN: Seen md5sum previously" >> "$DOCKERFILE"
  echo "" >> "$DOCKERFILE"

fi

done
rm "$tmp_seen"


echo -e '
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
' >> "$DOCKERFILE"


#exit 0
## Build using docker
$DOCKER build -t $USER/naive_atlas:${VERSION} .
$DOCKER run $USER/naive_atlas:${VERSION} naive_atlas --help

## Print commands to check the build
echo -e "

# Run these commands once you are happy with the containerize version
$DOCKER push $USER/naive_atlas:${VERSION}
docker image ls
docker image rm XXXX


# Once it is uploaded, test using singularity
singularity pull naive_atlas_${VERSION}.sif docker://$USER/naive_atlas:${VERSION}
singularity exec naive_atlas_${VERSION}.sif naive_atlas --help
rm naive_atlas_${VERSION}.sif

# OR build a SIF file directly from a local docker image
singularity build naive_atlas_${VERSION}.sif docker-daemon://$USER/naive_atlas:${VERSION}
singularity exec naive_atlas_${VERSION}.sif naive_atlas --help
rm naive_atlas_${VERSION}.sif


docker image ls
docker image rm XXXX
"


