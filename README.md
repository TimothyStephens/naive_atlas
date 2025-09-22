# Metagenome-Naive_Atlas

Metagenome-naive_atlas is a easy-to-use metagenomic pipeline based on snakemake. It handles all steps from QC, Assembly, Binning, to Annotation, and is designed to be as domain agnostic as possible (i.e., to assemble and bin eukaryotes, prokaryotes, and viruses).
naive_atlas is built upon the [ATLAS](https://github.com/metagenome-atlas/atlas) workflow, with enhancments from [VEBA](https://github.com/jolespin/veba) which allow it to identify MAGs from all domains.
All credit should go to the original authors of both workflows. 


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


# Citations

> ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.  
> Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A.   
> BMC Bioinformatics 21, 257 (2020).  
> doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)

> Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing.
> Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL.
> Nucleic Acids Res. 2024 Jun 22:gkae528.
> doi: [10.1093/nar/gkae528](https://doi.org/10.1093/nar/gkae528). PMID: 38909293.


# Sample file

The sample file used to tell the workflow about your samples has 4 required columns: "Reads_raw_R1", "Reads_raw_R2", "Assembler", "Bin_group"
 - `Reads_raw_R1`:    Path to first short read file (first mate of paired-end, single-end, or interleaved reads)
 - `Reads_raw_R2`:    Path to second short read file (second mate of paired-end, leave blank if single-end, or interleaved reads)
 - `Assembler`:       Assembler to us (see below)
 - `Bin_group`:       Groups to use when mapping read data for binning (all sample in a group will be mapped against each other, more samples helps binning, but massivly increases runtime)

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


Optional extra columns:
 - `Reads_raw_Long`                     Path to long reads (PacBio or Nanopore)
 - `Interleaved`                        Is the R1 short read file interleaved?                                      Options: `True` or `False`; False by default
 - `DeDuplicate`                        Should the sample's reads have depuplication run on itbefore use?           Options: `True` or `False`; True by default
 - `Quality_filter`                     Should the sample's reads have quality filtering run on it before use?      Options: `True` or `False`; True by default
 - `Remove_contaminants`                Should the sample's reads have contaminant sequences removed before use?    Options: `True` or `False`; True by default
 - `Normalize_reads_before_assembly`    Should the sample's reads be normalized before assembly?                    Options: `True` or `False`; True by default
 - `Error_correction_before_assembly`   Should the sample's reads be error corrected before assembly?               Options: `True` or `False`; True by default



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


# Known bugs

 - Parts of the workflow will randomly fail if you have packages (namely `numpy`) installed in `$USER/.local/lib/python*`. Python within the snakemake jobs will try and use the packages in theis directory over the ones installed in their conda environment, causing version issues with packages like `numpy`. The easiest way is to remove this whole directory to force python to always use it intended local conda env.





