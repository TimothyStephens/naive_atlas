# Metagenome-Naive_Atlas

Metagenome-naive_atlas is a easy-to-use metagenomic pipeline based on snakemake. It handles all steps from QC, Assembly, Binning, to Annotation, and is designed to be as domain agnostic as possible (i.e., to assemble and bin eukaryotes, prokaryotes, and viruses).
naive_atlas is built upon the [ATLAS](https://github.com/metagenome-atlas/atlas) workflow, with enhancements from [VEBA](https://github.com/jolespin/veba) which allow it to identify MAGs from all domains.
All credit should go to the original authors of both workflows. 





## Installation

You can start using naive_atlas with the following commands:
```
# Setup
git clone https://github.com/TimothyStephens/naive_atlas.git
cd naive_atlas

mamba env create -n naive_atlas --file naive_atlasenv.yml
conda activate naive_atlas

$CONDA_PREFIX/bin/pip3 install --prefix $CONDA_PREFIX --editable .
```

You can now use the workflow by just loading the conda env and running the `naive_atlas` command.
```bash
conda activate naive_atlas
naive_atlas --help
```





## Setup

To run `naive_atlas` you will need two files: a config file (`config.yaml`) and a samples manifest (`samples.tsv`).
Templates of these files can be created using the following command.
```bash
naive_atlas init --working-dir /path/to/your/project/dir
```
This will create a `config.yaml` and `samples.tsv` files in `/path/to/your/project/dir`.
You will need to manually edit `samples.tsv` with your own samples information.
The `config.yaml` will also need to be modified if you have reference host genome(s) that you would like to map against, or you wish to run additional annotation steps.
If you dont have reference host genome(s) and are happy with the minimal (fast) annotation steps, then this file can be left untouched.

### Sample file

The `samples.tsv` file is where you tell the workflow about your samples (type of reads, assembler to use, etc.). It will have 5 column: "Reads_raw_R1", "Reads_raw_R2", "Reads_raw_Long", "Assembler", "Bin_group"
 - `Reads_raw_R1`:    Path to first short read file (first mate of paired-end, single-end, or interleaved reads) (only required column)
 - `Reads_raw_R2`:    Path to second short read file (second mate of paired-end, leave blank if single-end, or interleaved reads)
 - `Reads_raw_Long`:  Path to long reads (PacBio or Nanopore)
 - `Assembler`:       Assembler to us (see below)
 - `Bin_group`:       Groups to use when mapping read data for binning (all sample in a group will be mapped against each other, more samples helps binning, but massively increases runtime)

The possible options for the `Assembler` column are:
 - `megahit`                Short reads (PE or SE)

 - `spades                  Short reads (PE or SE)
 - `spades-pacbio-raw`      Short reads (PE or SE) + PacBio regular CLR reads (<20% error)
 - `spades-pacbio-corr`     Short reads (PE or SE) + PacBio reads that were corrected with other methods (<3% error)
 - `spades-pacbio-hq`       Short reads (PE or SE) + PacBio HiFi reads (<1% error)
 - `spades-nanopore-raw`    Short reads (PE or SE) + ONT regular reads, pre-Guppy5 (<20% error)
 - `spades-nanopore-corr`   Short reads (PE or SE) + ONT reads that were corrected with other methods (<3% error)
 - `spades-nanopore-hq`     Short reads (PE or SE) + ONT high-quality reads (<1% error)

 - `flye-pacbio-raw`        Short reads (PE or SE) + PacBio regular CLR reads (<20% error)
 - `flye-pacbio-corr`       Short reads (PE or SE) + PacBio reads that were corrected with other methods (<3% error)
 - `flye-pacbio-hq`         Short reads (PE or SE) + PacBio HiFi reads (<1% error)
 - `flye-nanopore-raw`      Short reads (PE or SE) + ONT regular reads, pre-Guppy5 (<20% error)
 - `flye-nanopore-corr`     Short reads (PE or SE) + ONT reads that were corrected with other methods (<3% error)
 - `flye-nanopore-hq`       Short reads (PE or SE) + ONT high-quality reads (<1% error)

 - `metamdbg-pacbio-hq`     PacBio HiFi reads (<1% error)
 - `metamdbg-nanopore-hq`   ONT high-quality reads (<1% error)


Optional extra columns (normally these are not required and will be added by the workflow if not in the `samples.tsv` file):
 - `Interleaved`                        Is the R1 short read file interleaved?                                      Options: `True` or `False`; False by default
 - `DeDuplicate`                        Should the sample's reads have depuplication run on itbefore use?           Options: `True` or `False`; True by default
 - `Quality_filter`                     Should the sample's reads have quality filtering run on it before use?      Options: `True` or `False`; True by default
 - `Remove_contaminants`                Should the sample's reads have contaminant sequences removed before use?    Options: `True` or `False`; True by default
 - `Normalize_reads_before_assembly`    Should the sample's reads be normalized before assembly?                    Options: `True` or `False`; True by default
 - `Error_correction_before_assembly`   Should the sample's reads be error corrected before assembly?               Options: `True` or `False`; True by default


When you run the workflow for the first time a file called `naive_atlas.samples.run_info.tsv` will be created.
This is the full sample table (with all missing columns automatically added by the workflow) and the actual table that the workflow uses in the background when it runs.
Normally you dont have to worry about this file, although it can be helpful when assessing how your samples ran or why certain steps are being skipped or run for each sample.





### Config file

The `config.yaml` file contains parameters that are most useful for the user when running the workflow.
This file extends (and overwrites) the full config file in `naive_atlas/config/default_config.yaml`; the `default_config.yaml` file has a lot of parameters that the user normall does not need to worry about, so can be safely ignored.


#### Execution parameters

Runtime and memory parameters to use for different parts of the workflow.
This can mostly be ignored and is mainly used for execution on SLURM or other schedule systems.
Although the memory params are useful for local execution when using `--max-mem GBs`.
```yaml
########################
# Execution parameters
########################
#Runtime only for cluster execution

# BBMAPs tools threads and memory - needs relatively high amount of memory: normalize_reads, error_correction, merge_pairs
large_memory: 150 #GB
large_threads: 46
large_runtime: 24 #hr
# Assembly threads and memory: run_megahit, run_spades
assembly_memory: 350 #GB
assembly_threads: 60
assembly_runtime: 48 #hr
# All other major analysis rules Simple text parsing jobs
simplejob_memory: 20 #GB
simplejob_threads: 8
simplejob_runtime: 24 #hr
# Simple text parsing jobs
mem: 2 #GB
threads: 8
runtime: 4 #hr
```


#### Paths to reference genomes and databases

Paths to databases, tmp dir, and reference host datasets.
Generally, the `init` command will set the paths correctly, although the location of the `database_dir` can be changed if needed.
Also, the user should/can add their reference genome(s) to `contaminant_references` as required.
```yaml
########################
# Paths to reference genomes and databases
########################
# directory where databases are downloaded with 'atlas download'
database_dir:
  /project/databases

tmpdir:
  /project/tmp

# used to trim adapters from reads and read ends
preprocess_adapters:
  /project/databases/adapters.fa

#contamination references can be added such that -- key: /path/to/fasta
contaminant_references:
  PhiX:
    /project/databases/phiX174_virus.fa
  #Host_A:
  #  /project/host_refs/Host_A.fa
  #Host_B:
  #  /project/host_refs/Host_B.fa
```


#### Annotations

Types of annotation tools to run on predicted genes and genomes (MAGs and unbinned contigs).
The ones that are uncommented are the default minimal (fast) annotation approaches.
If the user wants to run additional tools, uncomment them as required.
```yaml
########################
# Annotations
########################
genome_annotations:
- gtdb_tree
- gtdb_taxonomy
#- dram
#- dram_unbinned
#- metaeuk
#- metaeuk_unbinned
#- mmseqs2_easy_taxonomy
#- mmseqs2_easy_taxonomy_unbinned

gene_annotations:
- eggNOG
#- eggNOG_unbinned
- dram
#- dram_unbinned
#- mmseqs2
#- mmseqs2_unbinned
```





## Run Workflow

Run the workflow on the user modified `config.yaml` and `samples.tsv` files.

To run just the `qc` stage.
This command will run as many jobs as it can until it hits 72 total threads or 400GB of RAM. 
The `--apptainer-args` parameter is required since some tools are run in containers and we need to define which directories it can see. You will need to change `/project` to whatever makes sense for your local system.
```bash
naive_atlas run qc \
  -c config.yaml \
  -j 72 --max-mem 400 \
  --apptainer-args '"--bind /project:/project"' --nt
```

To run all stages.
```bash
naive_atlas run all \
  -c config.yaml \
  -j 72 --max-mem 400 \
  --apptainer-args '"--bind /project:/project"' --nt
```





## Sample Download

`naive_atlas` does not have an inbuilt function to download datasets from NCBI.
HOWEVER, this probrlem has already been solved by [Kingfisher](https://wwood.github.io/kingfisher-download), which is automatically downloaded when setting up `naive_atlas`.
To download one or multiple SRA runs.
```bash
kingfisher get -r SRR14611058 SRR14611059 -m ena-ascp ena-ftp prefetch --output-directory SRA -f "fastq.gz"
```
To download all runs from one or multiple BioProjects.
```bash
kingfisher get -p PRJNA731596 PRJNA694677 -m ena-ascp ena-ftp prefetch --output-directory SRA -f "fastq.gz"
```
These commands will download the resulting `fastq.gz` files into the directory `SRA/`.
Once you have downloaded these files, they can be added to your `samples.tsv` file.





## Cluster Execution

To run `naive_atlas` on a SLURM system, you can use the cluster profile in `naive_atlas/config/slurm_profile`.
Within this directory is a submission script `run_naive_atlas.sh`, which is a template for the job which will run `naive_atlas` (the job from which all other jobs will be submitted and managed).
Since shared SLURM systems often have preemption enabled or maximum walltimes, which will result in the `naive_atlas` job from being canceled before the workflow can finish, this script is designed to resubmit its self if preempted or killed becuase of walltime.





## Output files

Output files will be in `genomes/`.
Not all annotation files will be present if they are not selected to in the `config.yaml` file.
The sub result directories:

- `genomes/alignments/` Alignment files.
  - `all_contigs.fa` All MAG scaffolds used for mapping.
  - `bams/` Alignment BAM files against `all_contigs.fa`, one per sample.

- `genomes/annotations/` Genome annotation files.
  - `genomes/annotations/genomes/` Annotation files for each MAG.
    - `dram/annotations.tsv` DRAM annotations.
    - `dram/kegg_modules.tsv` DRAM KEGG annotations.
    - `mmseqs2/` MMSEQS2 annotation of genes from each MAG.
      - `*.easy_taxonomy_result_report` Kraken-like report from MMSEQS2 easy-taxonomy, can be visualized using [pavian](https://github.com/fbreitwieser/pavian).
    - `genes/eggNOG.tsv.gz` EggNOG-mapper results.
    - `genes/dram` DRAM results.
    - `metaeuk/` MetaEuk taxonomic annotations.
    - `metaeuk_contig_predictions.tsv` Combined per contig MetaEuk results (from all MAGs).
    - `metaeuk_mag_predictions.tsv` Combined per MAG MetaEuk results.
    - `taxonomy/gtdb_taxonomy.tsv` GTDB Taxonomic results.
    - `tree/` Phylogeny results.
  - `genomes/annotations/unbinned/` Annotation files for each MAG.
    - `dram/annotations.tsv` DRAM annotations.
    - `dram/kegg_modules.tsv` DRAM KEGG annotations.
    - `genes/mmseqs2/` MMSEQS2 annotation of genes from each sample.
      - `*.easy_taxonomy_result_report` Kraken-like report from MMSEQS2 easy-taxonomy, can be visualized using [pavian](https://github.com/fbreitwieser/pavian).
    - `genes/eggNOG.tsv.gz` EggNOG-mapper results.
    - `genes/dram` DRAM results.
    - `metaeuk/` MetaEuk taxonomic annotations.
    - `metaeuk_contig_predictions.tsv` Combined per contig MetaEuk results (from all samples unbinned scaffolds).
    - `metaeuk_mag_predictions.tsv` Combined per sample MetaEuk results.

- `genomes/clustering/` Results from MAG clustering, can be used to trace the MAGs from each sample back to the final MAGs in `genomes/genomes/`.

- `genomes/coverage/` MAG coverage results. NOTE: Only the full QC reads are used for mapping, so reads filtered becuase they aligned to the supplied "contaminant" datasets will not be present in the results. If you want a host vs. MAG coverage analysis you will have to run it yourself.
  - `coverage.tsv.gz` Coverage results from `coverm`. This is the main results file for plotting.
  - `read_stats.tsv` Mapping stats from coverage analysis. Can be used to assess how many reads mapped against your MAGs.

- `genomes/genes/` Predicted coding and non-coding genes.
  - `genomes/genes/genomes` Predicted genes in each MAG. Some files are specific to prokaryotes vs. eukaryotes.
    - `*.faa` Predicted gene protein sequences.
    - `*.fna` Predicted gene nucleotide (CDS) sequences.
    - `*.gff3` Predicted gene feature coordinates along scaffolds in MAG.
    - `*.ncRNA.fna` Predicted non-coding genes nucleotide (CDS) sequences.
    - `*.ncRNA.gff3` Predicted non-coding genes coordinates along scaffolds in MAG.
    - `*.other.fna` Predicted other features (e.g., crispr repeats) nucleotide (CDS) sequences.
    - `*.other.gff3` Predicted other features (e.g., crispr repeats) coordinates along scaffolds in MAG.
    - `*.rRNA.fna` Predicted rRNA genes nucleotide (CDS) sequences.
    - `*.rRNA.gff3` Predicted rRNA genes coordinates along scaffolds in MAG.
    - `*.tRNA.fna` Predicted tRNA genes nucleotide (CDS) sequences.
    - `*.tRNA.gff3` Predicted tRNA genes coordinates along scaffolds in MAG.
    - `*.tsv` Predicted gene information table.
    - `*.gene_stats.tsv` Gene stats for each group of MAGs (eukaryotes, prokaryotes, viruses, plasmids)

  - `genomes/genes/unbinned/` Predicted genes in each sample's unbinned scaffolds.
    - `*.faa` Predicted gene protein sequences.
    - `*.fna` Predicted gene nucleotide (CDS) sequences.
    - `*.gff3` Predicted gene feature coordinates along unbinned scaffolds.

- `genomes/genomes/` MAG scaffolds.
  - `*.fa` MAGs.
  - `*.genome_quality.tsv` MAG quality metrics.

- `genomes/unbinned/` Unbinned scaffolds from each sample.





## Known bugs

 - Parts of the workflow will randomly fail if you have packages (namely `numpy`) installed in `$USER/.local/lib/python*`. Python within the snakemake jobs will try and use the packages in theis directory over the ones installed in their conda environment, causing version issues with packages like `numpy`. The easiest way is to remove this whole directory to force python to always use it intended local conda env.





## Citations

If you use this workflow, please cite this GitHub repository and the following manuscripts to show support for the tools that this workflow is built on.

> ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.
> Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A.
> BMC Bioinformatics 21, 257 (2020).
> doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)

> Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing.
> Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL.
> Nucleic Acids Res. 2024 Jun 22:gkae528.
> doi: [10.1093/nar/gkae528](https://doi.org/10.1093/nar/gkae528). PMID: 38909293.


