# Metagenome-Naive_Atlas

[![Anaconda-Server Badge](https://anaconda.org/bioconda/metagenome-atlas/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/metagenome-atlas)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda )](https://anaconda.org/bioconda/metagenome-atlas)
[![Documentation Status](https://readthedocs.org/projects/metagenome-atlas/badge/?version=latest)](https://metagenome-atlas.readthedocs.io/en/latest/?badge=latest)
![Mastodon Follow](https://img.shields.io/mastodon/follow/109273833677404282?domain=https%3A%2F%2Fmstdn.science&style=social)
<!--[![follow on twitter](https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow)](https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd) -->


Metagenome-naive_atlas is a easy-to-use metagenomic pipeline based on snakemake. It handles all steps from QC, Assembly, Binning, to Annotation, and is designed to be as domain agnostic as possible (i.e., to assemble and bin eukaryotes, prokaryotes, and viruses).
All credit for this amazing workflow should go to the original authors. naive_atlas is simply a reconfiguration to remove filtering steps that select for only prokaryotes, allowing for eukaryote and viral bins to also be produced by the workflow.

![scheme of workflow](resources/images/atlas_list.png?raw=true)

You can start using naive_atlas with the following commands:
```
# Setup
git clone https://github.com/TimothyStephens/naive_atlas.git
cd naive_atlas

mamba env create -n naive_atlas --file naive_atlasenv.yml
conda activate naive_atlas

$CONDA_PREFIX/bin/pip3 install --prefix $CONDA_PREFIX --editable .

# Analysis
naive_atlas init --db-dir databases path/to/fastq/files
naive_atlas run all
```
naive_atlas does not have its own dedicted documentation, however, the atlas documentation is still highly relevent.

# Webpage

[metagenome-atlas.github.io](https://metagenome-atlas.github.io/)

# Documentation

https://metagenome-atlas.readthedocs.io/

[Tutorial](https://github.com/metagenome-atlas/Tutorial)

# Citation

> ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.  
> Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A.   
> BMC Bioinformatics 21, 257 (2020).  
> doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)


# Developpment/Extensions

Here are some ideas I work or want to work on when I have time. If you want to contribute or have some ideas let me know via a feature request issue.

- Optimized MAG recovery (e.g. [Spacegraphcats](https://github.com/spacegraphcats/spacegraphcats))
- Integration of viruses/plasmid that live for now as [extensions](https://github.com/metagenome-atlas/virome_atlas)
- Add statistics and visualisations as in [atlas_analyze](https://github.com/metagenome-atlas/atlas_analyze)
- Implementation of most rules as snakemake wrapper
- Cloud execution
- Update to new Snakemake version and use cool reports.

# Output files

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


