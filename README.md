# Metagenome-Naive_Atlas

Metagenome-naive_atlas is an easy-to-use metagenomic pipeline built on Snakemake. It handles all steps from QC, assembly, binning, gene prediction, functional annotation, and genome and gene quantification. It is designed to be domain agnostic, assembling and binning prokaryotes, (micro)eukaryotes, viruses, and plasmids.

naive_atlas is built upon the [ATLAS](https://github.com/metagenome-atlas/atlas) workflow, with enhancements from [VEBA](https://github.com/jolespin/veba) that enable it to recover MAGs from all domains of life.


## Setup

### Local Install

```bash
git clone https://github.com/TimothyStephens/naive_atlas.git
cd naive_atlas

mamba env create -n naive_atlas --file naive_atlasenv.yml
conda activate naive_atlas

$CONDA_PREFIX/bin/pip install --prefix $CONDA_PREFIX --editable .
```

### Cluster Install

For cluster execution (SLURM or generic), use the cluster environment file which includes executor plugins:

```bash
git clone https://github.com/TimothyStephens/naive_atlas.git
cd naive_atlas

mamba env create -n naive_atlas --file naive_atlasenv_cluster.yml
conda activate naive_atlas

$CONDA_PREFIX/bin/pip install --prefix $CONDA_PREFIX --editable .
```

### Docker

The `docker_build` script generates Dockerfiles and builds Docker images from the workflow's conda environments:

```bash
./docker_build
```

Use `--build_only` to generate the Dockerfiles without building. The script produces two images: `timothystephens/naive_atlas-envs` (conda environments) and `timothystephens/naive_atlas` (full workflow).


## Quick Start

```bash
naive_atlas init --db-dir databases path/to/fastq/files
naive_atlas run all
```

Run `naive_atlas run --help` for options.


## Pipeline Stages

The pipeline follows this DAG:

```
qc -> assembly -> binning -> genomes -> quantify_genomes -> strains
                              +-> genome_annotation
                                            +-> gene_prediction -> quantify_genes -> gene_annotation
```

Plus two independent stages:
- `screen` — rapid screening via skani (runs after QC)
- `download` — pre-download all reference databases

| Stage | Description |
|---|---|
| `download` | Pre-download all reference databases (~1.2 TB) |
| `qc` | Quality control: dedup, quality filter, contaminant removal, read stats |
| `assembly` | Assembly with megahit/spades/flye/metamdbg, contig filtering |
| `binning` | Multi-domain binning: prokaryotic (MetaBAT2, MaxBin2, MDMcleaner, CheckM2), eukaryotic (MetaBAT2, SemiBin, BUSco), viral (CheckV, GUNC) |
| `genomes` | Species-level clustering via skANI, CD-HIT deduplication |
| `quantify_genomes` | MAG coverage quantification (coverm) |
| `genome_annotation` | GTDB-Tk taxonomy, MetaEuk, MMseqs2 taxonomy |
| `gene_prediction` | Prokaryotic (Prodigal), eukaryotic (MetaEuk, MicroEuk), ncRNA (Barrnap, Infernal) |
| `quantify_genes` | Gene-level coverage quantification |
| `gene_annotation` | eggNOG-mapper, MMseqs2 easy-search |
| `strains` | Strain-level analysis via inStrain |
| `screen` | Rapid screening via skani against references (runs after QC) |
| `all` | Run complete pipeline |


## Citations

> ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.  
> Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A.   
> BMC Bioinformatics 21, 257 (2020).  
> doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)

> Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing.  
> Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL.  
> Nucleic Acids Res. 2024 Jun 22:gkae528.  
> doi: [10.1093/nar/gkae528](https://doi.org/10.1093/nar/gkae528). PMID: 38909293.


## Sample File (`samples.tsv`)

The sample file tells the workflow about your samples. It has 4 required columns: `Reads_R1`, `Reads_R2`, `Assembler`, `Bin_group`.

- `Reads_R1` — Path to first short read file (first mate of paired-end, single-end, or interleaved reads)
- `Reads_R2` — Path to second short read file (second mate of paired-end; leave blank if single-end or interleaved)
- `Assembler` — Assembler to use (see below)
- `Bin_group` — Group name for co-binning. All samples in a group are mapped against each other; more samples per group improves binning but significantly increases runtime.

### Assembler Options

- `megahit` — Short reads (PE or SE)
- `spades` — Short reads (PE or SE)
- `spades-pacbio-raw` — Short reads (PE only) + PacBio CLR reads (<20% error)
- `spades-pacbio-corr` — Short reads (PE only) + PacBio corrected reads (<3% error)
- `spades-pacbio-hq` — Short reads (PE only) + PacBio HiFi reads (<1% error)
- `spades-nanopore-raw` — Short reads (PE only) + ONT regular reads, pre-Guppy5 (<20% error)
- `spades-nanopore-corr` — Short reads (PE only) + ONT corrected reads (<3% error)
- `spades-nanopore-hq` — Short reads (PE only) + ONT high-quality reads (<1% error)
- `flye-pacbio-raw` — PacBio CLR reads (<20% error)
- `flye-pacbio-corr` — PacBio corrected reads (<3% error)
- `flye-pacbio-hq` — PacBio HiFi reads (<1% error)
- `flye-nanopore-raw` — ONT regular reads, pre-Guppy5 (<20% error)
- `flye-nanopore-corr` — ONT corrected reads (<3% error)
- `flye-nanopore-hq` — ONT high-quality reads (<1% error)
- `metamdbg-pacbio-hq` — PacBio HiFi reads (<1% error)
- `metamdbg-nanopore-hq` — ONT high-quality reads (<1% error)

### Optional Columns

- `Reads_Long` — Path to long reads (PacBio or Nanopore)
- `Interleaved` — Is the R1 short read file interleaved? Options: `True` or `False`; `False` by default
- `DeDuplicate` — Should deduplication be run on this sample's reads? Options: `True` or `False`; `True` by default
- `Quality_filter` — Should quality filtering be run on this sample's reads? Options: `True` or `False`; `True` by default
- `Remove_contaminants` — Should contaminant sequences be removed from this sample? Options: `True` or `False`; `True` by default
- `Normalize_reads_before_assembly` — Should reads be normalized before assembly? Options: `True` or `False`; `True` by default
- `Error_correction_before_assembly` — Should error correction be run before assembly? Options: `True` or `False`; `True` by default

# Configuration

## How Configuration Works

Running `naive_atlas init` generates a `config.yaml` in the project directory, copied from the package's `config/template_config.yaml`. This is the file users should edit.

The template exposes a small set of commonly modified parameters. The full configuration lives in `config/default_config.yaml`, which Snakemake merges with the user's `config.yaml` — **any key from the default config can be overridden** by adding it to your `config.yaml`.

You can also supply `--config-file path/to/other.yaml` on the CLI to use a different config file.

---

## Parameters (Default Sser `config.yaml`)

These are the parameters presented in `config.yaml` after running `naive_atlas init`. They cover the most common customization points without requiring knowledge of the full config system.

### `database_dir`

Path where reference databases are stored. All databases combined require approximately 1.2 TB of disk space.

* **Default:** `"/user/project/dir/databases"`
* **When to change:** Always — set this to the actual path where databases will be downloaded or already reside.

### `contaminant_references`

A dictionary mapping contaminant reference names to FASTA file paths. During QC, reads that match any of these references are removed. This is the primary mechanism for host-read decontamination (e.g., removing human reads from stool samples).

* **Default:** `{}` (no contaminants removed)
* **Format:** YAML dictionary; keys are arbitrary names, values are paths to FASTA files.
* **Commented out in the template** — uncomment and fill in to activate.

```yaml
contaminant_references:
  human: /path/to/human/genome.fa
  vector: /path/to/vector/sequences.fa
```

* **When to change:** Whenever your samples may contain reads from a host or other contaminant organism. Without this, all reads pass through QC (including host).

### `contaminant_references_include_phiX`

When `true`, the PhiX control virus sequence is automatically included as a contaminant reference. This is added in addition to any user-specified `contaminant_references`. Controlled at the Snakemake level via `resources.smk`, not in the template file itself.

* **Default:** `true`
* **When to change:** Set to `false` only if you deliberately want to retain PhiX reads (rare).

### `veba_prokaryotic.mdmcleaner_skip`

Controls whether the MDMcleaner dereplication step runs during prokaryotic binning. MDMcleaner removes chimeric and contaminant contigs from bins using GTDB reference data, but can be extremely slow on rich prokaryotic microbiomes (e.g., soil, activated sludge).

* **Default:** `false` (MDMcleaner runs)
* **When to change:** Set to `true` to skip MDMcleaner and significantly reduce runtime. The trade-off is that bins may contain more chimeric or contaminant contigs.

### `veba_prokaryotic.mdmcleaner_options`

Extra command-line options passed directly to the MDMcleaner invocation.

* **Default:** `''` (no extra options)
* **Example:** `'--fast_run'` to reduce MDMcleaner runtime at the cost of slightly less thorough cleaning.

### `genome_annotations`

A list of genome annotation tasks to run on binned MAGs. Each entry corresponds to a Snakemake rule. The template enables GTDB taxonomy by default and comments out the rest.

* **Default in template:**
  ```yaml
  genome_annotations:
  - gtdb_tree
  - gtdb_taxonomy
  ```
* **Available options:**
  - `gtdb_tree` — Build a phylogenetic tree with GTDB-Tk (lightweight).
  - `gtdb_taxonomy` — Classify MAGs with GTDB-Tk taxonomy (lightweight).
  - `metaeuk` — MetaEuk gene prediction + taxonomic classification on binned MAGs (resoruce intensive).
  - `metaeuk_unbinned` — MetaEuk on unbinned scaffolds (resource intensive).
  - `mmseqs2_easy_taxonomy` — Kraken-like taxonomy assignment via MMseqs2 easy-taxonomy on binned MAGs (results viewable in [Pavian](https://github.com/fbreitwieser/pavian)) (resource intensive).
  - `mmseqs2_easy_taxonomy_unbinned` — MMseqs2 easy-taxonomy on unbinned scaffolds (resource intensive).
* **When to change:** Uncomment additional entries to enable more annotation layers. The `gtdb_taxonomy` + `mmseqs2_easy_taxonomy` combination gives broad taxonomic coverage.

### `gene_annotations`

A list of gene functional annotation tasks to run on predicted genes from MAGs.

* **Default in template:**
  ```yaml
  gene_annotations:
  - eggNOG
  - eggNOG_unbinned
  ```
* **Available options:**
  - `eggNOG` — eggNOG-mapper functional annotation on binned MAG genes (resource intensive).
  - `eggNOG_unbinned` — eggNOG-mapper on unbinned scaffold genes (resource intensive).
  - `mmseqs2_easy_search` — MMseqs2 easy-search against the configured `mmseqs2_database` on binned MAG genes (resource intensive).
  - `mmseqs2_easy_search_unbinned` — MMseqs2 easy-search on unbinned scaffold genes (resource intensive).
* **When to change:** Uncomment `mmseqs2_easy_search` entries to get additional functional annotations beyond eggNOG. The choice of database (`mmseqs2_database`) controls coverage breadth — see below.

---

## Advanced Parameters (Default Config Overrides)

The following parameters are **not** in the template file but can be added to your `config.yaml` to override the defaults from `config/default_config.yaml`. These are useful when the default behavior does not fit your data or cluster environment.

### QC Tuning

| Parameter | Default | Description |
|---|---|---|
| `preprocess_qtrim` | `"rl"` | Quality trimming mode (`"rl"` = trim both ends). |
| `preprocess_adapters` | `"workflow/data/adapters.fa"` | Path to a custom adapter FASTA file. Override if your library uses non-standard adapters. |
| `minimum_contig_length` | `1000` | Minimum contig length (bp) to retain after assembly. Increase to filter out more assembly fragments. |
| `preprocess_minimum_passing_read_length` | `51` | Minimum read length after quality trimming. Reads shorter than this are discarded. |
| `preprocess_minimum_base_quality` | `10` | Minimum Phred quality score for base trimming. |

> **Note:** Read deduplication and read normalization before assembly are controlled via the `DeDuplicate` and `Normalize_reads_before_assembly` columns in your `samples.tsv` file, not via `config.yaml`.

### Assembly Tuning

| Parameter | Default | Description |
|---|---|---|
| `minimum_average_coverage` | `1` | Minimum average coverage for contigs to be retained. Increase to require stronger read support. |
| `megahit_min_count` | `2` | Minimum k-mer count for MeGAHit. `2` for metagenomes, `3` for high-coverage genomes. |
| `megahit_k_max` | `121` | Maximum k-mer size for MeGAHit. Larger values improve assembly but use more memory. |
| `error_correction_kmer` | `31` | K-mer size for Tadpole error correction. Larger values (e.g., `62`) improve accuracy but use more memory. |

### Binning Thresholds

Each VEBA domain has independent quality thresholds. Override these to accept lower-quality bins (higher recall, lower precision) or stricter bins (higher precision, lower recall).

| Parameter | Default | Description |
|---|---|---|
| `veba_prokaryotic.minimum_contig_length` | `1500` | Minimum contig length for prokaryotic binning. |
| `veba_prokaryotic.minimum_genome_length` | `150000` | Minimum prokaryotic MAG length (bp). |
| `veba_prokaryotic.checkm2_completeness` | `50` | Minimum CheckM2 completeness (%) for a prokaryotic MAG. |
| `veba_prokaryotic.checkm2_contamination` | `10` | Maximum CheckM2 contamination (%) for a prokaryotic MAG. |
| `veba_eukaryotic.minimum_genome_length` | `2000000` | Minimum eukaryotic MAG length (bp). |
| `veba_eukaryotic.busco_completeness` | `10` | Minimum BUSco completeness (%) for a eukaryotic MAG. |
| `veba_eukaryotic.busco_contamination` | `10` | Maximum BUSco contamination (%) for a eukaryotic MAG. |
| `veba_viral.minimum_genome_length` | `2500` | Minimum viral/plasmid genome length (bp). |
| `veba_viral.minimum_score` | `0.7` | Minimum VIBRANT score to flag a sequence as viral. |

### Genome Dereplication

| Parameter | Default | Description |
|---|---|---|
| `genome_dereplication.ANI` | `0.95` | ANI threshold for genome clustering. `0.95` ≈ species-level; lower values cluster more broadly. |
| `genome_dereplication.overlap` | `0.2` | Minimum fraction of the shorter genome that must align for clustering. |

### Quantification

| Parameter | Default | Description |
|---|---|---|
| `coverm_params` | `'--min-read-percent-identity 95 --min-covered-fraction 0 --min-read-aligned-percent 70'` | Extra parameters passed to coverm for genome coverage. |
| `coverm_stats` | `'relative_abundance mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm'` | Coverage statistics computed by coverm. |

### MMseqs2 Database

| Parameter | Default | Description |
|---|---|---|
| `mmseqs2_database` | `"UniRef100"` | Database for MMseqs2 taxonomy and functional searches. Options: `UniRef100`, `UniRef90`, `UniRef50`, `UniProtKB`, `UniProtKB/TrEMBL`, `UniProtKB/Swiss-Prot`, `NR`, `GTDB`. Lower UniRef clustering (e.g., `UniRef50`) is faster but less specific; `NR` is broadest but slowest. |

### Resource Constraints

Override these to match your cluster's node limits or to restrict resource usage.

| Parameter | Default | Description |
|---|---|---|
| `max_mem_mb` | `1500000` | Absolute maximum memory (MB) any single job may request. Lower this if your cluster has smaller maximum node memory. |
| `simple_job_mem_mb` | `5000` | Memory (MB) for simple jobs (report generation, stats, downloads). |
| `simple_job_threads` | `8` | CPU cores for simple jobs. |
| `simple_job_time_min` | `1440` | Walltime (minutes) for simple jobs. |

### SLURM Partition Routing

The `partition_ranges` key maps memory thresholds to SLURM partition names. Jobs are routed to the first partition whose threshold their memory request does not exceed.

```yaml
partition_ranges:
  main-redhat: 160000    # Jobs ≤ 160 GB go here
  mem-redhat: 1500000    # Jobs > 160 GB go here
```

Override this to match your cluster's actual partition names and memory tiers. If a job's memory request exceeds all thresholds, it falls through to the last partition.

---

### Dynamic Resource Scaling

Every Snakemake rule in naive_atlas specifies a `resources` block with `base_mem_mb`, `scaling_factor`, and `threads`. Memory is computed dynamically from input file sizes:

```
mem_mb = min(base_mem_mb + (input_size_mb × scaling_factor), max_mem_mb)
```

Example per-rule specifications:

| Rule | Threads | Base Memory (MB) | Scaling Factor |
|---|---|---|---|
| Deduplicate reads | 24 | 8000 | 3.5 |
| Run assembly | 24 | 12000 | 6.0 |
| Quality filter | 24 | 8000 | 2.0 |
| Decontamination | 24 | 12000 | 1.5 |
| Binning | 24 | 40000 | 2.0 |
| Gene prediction | 24 | 80000 | 2.0 |

A 1 GB input file to the assembly rule requests approximately 12000 + (1024 × 6) = ~18144 MB. To override per-rule resources, add entries under the `resources` key in your `config.yaml` following the same structure as `config/default_config.yaml`. If a rule fails, then it will be rerun with additional memory since this is the most common reason for a job failing (especially on a cluster system). Rerun rules memory will increase by X, where X is the number of attempt: e.g., first attemp: 20Gb, sencond attempt: 40GB, etc.


## Cluster Execution

To run on a cluster, use `--cluster-type` with the `run` command:

```bash
naive_atlas run all --cluster-type slurm
naive_atlas run all --cluster-type generic
```

Key cluster options:

- **`--cluster-type`** — Type of cluster: `slurm` or `generic`.
- **`--local-cores`** — Number of local cores per job. Defaults to max system cores.
- **`--jobs`** — Maximum number of parallel jobs. Default: `250`.
- **`--cores`** — Number of cores (used when `--cluster-type` is NOT set). Defaults to max system cores.
- **`--profile`** — Snakemake profile for cluster execution.
- **`--cluster-slurm-params`** — Parameters for SLURM cluster execution.
- **`--sdm`** — Software deployment method: `apptainer` or `conda`.


## Output Files

Output files are written to the working directory. Not all annotation files will be present if they were not selected in the `config.yaml`.

- **`genomes/genomes/`** — Representative MAGs (FASTA + quality metrics)
- **`genomes/genes/`** — Predicted genes (from MAGs/sequences in `genomes/` and `unbinned/`)
- **`genomes/coverage/`** — MAG coverage (coverm results)
- **`genomes/annotations/`** — Taxonomy, functional annotations
- **`genomes/alignments/`** — BAM alignments
- **`genomes/unbinned/`** — Unbinned scaffolds
- **`genomes/clustering/`** — MAG clustering results
- **`genomes/strains/`** — inStrain results
- **`genomes/screen/`** — skani screening results
- **`reports/`** — HTML reports (QC, assembly, binning, coverage)
- **`samples/`** — Per-sample intermediate files (QC, assembly, binning)
- **`stats/`** — Read statistics
- **`logs/`** — Job logs
- **`tmp/`** — Temporary files

### Coverage Note

MAG coverage results use only the full QC reads, so reads filtered because they aligned to contaminant references will not be present. For host vs. MAG coverage analysis, you will need to run it yourself.


## Known Issues

- Parts of the workflow may fail if you have packages (namely `numpy`) installed in `$USER/.local/lib/python*`. Python within Snakemake jobs will try to use the packages in that directory over the ones installed in their conda environment, causing version issues. The easiest fix is to remove that directory to force Python to always use its local conda environment.
