####
#### Script to automatically containerize and build naiveATLAS into docker
####
set -euo pipefail


## Env vars
USER="timothystephens"
DOCKER=$(which docker)
DOCKERFILE="dockerfile"
VERSION=$(git describe --tags --dirty --always)


## Build dockerfile
echo -e "# Dockerfile to containerize naiveATLAS workflow
FROM condaforge/mambaforge:latest


## (1/7) Set environment variables
# So we dont have to interactivly configure tzdata
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
ENV DEBIAN_FRONTEND=noninteractive

ENV CONDA_DIR=/conda-envs
ENV PATH=/conda-envs/naive_atlas/bin:\$PATH


RUN apt-get update && \\
    apt-get install -y \\
      zlib1g zlib1g-dev build-essential gcc && \\
    rm -rf /var/lib/apt/lists/*


## (2/7) Install Singulairty
RUN apt-get update && \\
    apt-get install -y --no-install-recommends \\
      software-properties-common \\
      gnupg \\
      ca-certificates \\
    && add-apt-repository -y ppa:apptainer/ppa \\
    && apt-get update \\
    && apt-get install -y apptainer \\
    && rm -rf /var/lib/apt/lists/*


## (3/7) Create and set workflow working directory
WORKDIR /app
COPY . .
RUN mamba env create --prefix /conda-envs/naive_atlas --file naive_atlasenv.yml
SHELL [\"conda\", \"run\", \"-p\", \"/conda-envs/naive_atlas\", \"/bin/bash\", \"-c\"]
RUN /conda-envs/naive_atlas/bin/pip3 install --prefix /conda-envs/naive_atlas --editable .
SHELL [\"/bin/sh\", \"-c\"]


## (4/7) Install each workflow package" > "$DOCKERFILE"

RUN=""
tmp_seen=$(mktemp)
for YAML in workflow/envs/*.yaml;
do

MD5SUM=$(md5sum "$YAML" | awk '{print $1}')

# Check if we have seen this md5sum previously
if ! grep -qFx "$MD5SUM" "$tmp_seen"; then
  # Not seen

echo -e "# Conda environment:
#   source: $YAML
#   prefix: /conda-envs/$MD5SUM" >> "$DOCKERFILE"
awk '{print "#   "$0}' "$YAML"   >> "$DOCKERFILE"
echo -e "COPY $YAML /conda-envs/$MD5SUM.yaml" >> "$DOCKERFILE"
echo "" >> "$DOCKERFILE"

# Need to run mamba create as each RUN creates a new layer, which is a lot of diskspace overhead and will cause issues if done separatly.
NL=$'\n'
if [ -z "$RUN" ]; then
  RUN="RUN mamba env create --prefix /conda-envs/$MD5SUM --file /conda-envs/$MD5SUM.yaml && \\"
else
  RUN="${RUN}${NL}    mamba env create --prefix /conda-envs/$MD5SUM --file /conda-envs/$MD5SUM.yaml && \\"
fi

  echo "$MD5SUM" >> "$tmp_seen"

else

echo -e "# Conda environment:
#   source: $YAML
#   prefix: /conda-envs/$MD5SUM" >> "$DOCKERFILE"
awk '{print "#   "$0}' "$YAML"   >> "$DOCKERFILE"
echo -e "# NOTRUN: Seen md5sum previously" >> "$DOCKERFILE"
echo "" >> "$DOCKERFILE"

fi

done
RUN="${RUN}${NL}    mamba clean -afy"
rm "$tmp_seen"
echo -e "$RUN" >> "$DOCKERFILE"



echo -e '

## (5/7) Download Singulairty/Apptainer SIF files
RUN mkdir -p /containers-sif' >> "$DOCKERFILE"

RUN=""
tmp_seen=$(mktemp)
for URI in $(cat workflow/rules/*.smk \
  | grep 'docker:' \
  | awk '{print $1}' \
  | sed -e 's/,//'); do

HASH=$(echo -n "$URI" | sha256sum | awk '{print $1}')
SIF="/containers-sif/${HASH}.sif"

if ! grep -qFx "$HASH" "$tmp_seen"; then
  # Not seen
  NL=$'\n'
  if [ -z "$RUN" ]; then
    RUN="RUN apptainer pull --name $SIF $URI && \\"
  else
    RUN="${RUN}${NL}    apptainer pull --name $SIF $URI && \\"
  fi
  
  echo "$HASH" >> "$tmp_seen"
fi

done

RUN="${RUN}${NL}    apptainer cache clean --force"
rm "$tmp_seen"
echo -e "\n$RUN" >> "$DOCKERFILE"


echo -e '

## (6/7) Verify the installation by checking the version
RUN naive_atlas --help


## (7/7) Set the default command to show MaxBin2 help
CMD ["naive_atlas", "--help"]
' >> "$DOCKERFILE"


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
singularity pull naive_atlas_v${VERSION}.sif docker://$USER/naive_atlas:${VERSION}
singularity exec naive_atlas_v${VERSION}.sif naive_atlas --help
rm naive_atlas_v${VERSION}.sif

docker image ls
docker image rm XXXX
"


