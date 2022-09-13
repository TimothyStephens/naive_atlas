## dRep
localrules:
    get_all_bins,
    all_contigs2bins,


rule get_all_bins:
    input:
        bins=expand(
            "{sample}/binning/{binner}/bins",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        temp(directory("genomes/all_bins")),
    run:
        os.mkdir(output[0])
        from utils.io import symlink_relative

        for bin_folder in input.bins:

            fasta_files = [f for f in os.listdir(bin_folder) if f.endswith(".fasta")]
            symlink_relative(fasta_files, bin_folder, output[0])


rule all_contigs2bins:
    input:
        expand(
            "{sample}/binning/{binner}/cluster_attribution.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        temp("genomes/clustering/all_contigs2bins.tsv.gz"),
    run:
        from utils.io import cat_files

        cat_files(input, output[0], gzip=True)


localrules:
    merge_checkm_for_dereplication,
    merge_checkm,


rule merge_checkm_for_dereplication:
    input:
        "reports/genomic_bins_{binner}.tsv".format(binner=config["final_binner"]),
    output:
        temp("genomes/quality.csv"),
    run:
        import pandas as pd

        D = pd.read_csv(input[0], index_col=0, sep="\t")


        D.index.name = "genome"
        
        D.columns = D.columns.str.lower()
        D = D[["completeness", "contamination"]]
        D.to_csv(output[0])


rule merge_checkm:
    input:
        completeness=expand(
            "{sample}/binning/{binner}/checkm/completeness.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
        taxonomy=expand(
            "{sample}/binning/{binner}/checkm/taxonomy.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
        markers=expand(
            "{sample}/binning/{binner}/checkm/storage/tree/concatenated.fasta",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        checkm="genomes/checkm/checkm_all_bins.tsv",
        markers="genomes/checkm/all_bins_markers.fasta",
    run:
        import pandas as pd
        import shutil

        D = []

        for i in range(len(SAMPLES)):
            df = pd.read_csv(input.completeness[i], index_col=0, sep="\t")
            df = df.join(pd.read_csv(input.taxonomy[i], index_col=0, sep="\t"))
            D.append(df)

        D = pd.concat(D, axis=0)
        D.to_csv(output.checkm, sep="\t")

        with open(output.markers, "wb") as fout:
            for fasta in input.markers:
                shutil.copyfileobj(open(fasta, "rb"), fout)


def get_dereplication_arguments(key):
    "key = completeness or contamination"

    if "filter" in config["genome_dereplication"]:
        logger.info(
            "The configuration for dereplication has changed. For now I will try to be backward compatible, but please update your config file. https://metagenome-atlas.readthedocs.io/en/latest/usage/configuration.html#detailed-configuration"
        )

        value = config["genome_dereplication"]["filter"][key]

    else:
        if key == "completeness":
            return config["genome_min_completeness"]
        elif key == "contamination":
            return config["genome_max_contamination"]
        else:
            raise ValueError("key must be completeness, contamination")


rule dereplication:
    input:
        dir="genomes/all_bins",
        quality="genomes/quality.csv",
    output:
        dir=temp(directory("genomes/dereplicated_genomes")),
        mapping_file=temp("genomes/clustering/allbins2genome_oldname.tsv"),
    threads: config["threads"]
    log:
        "logs/genomes/dereplication.log",
    conda:
        "../envs/galah.yaml"
    params:
        ANI=config["genome_dereplication"]["ANI"],
        quality_formula=config["genome_dereplication"]["quality_formula"],
        min_completeness=get_dereplication_arguments("completeness"),
        max_contamination=get_dereplication_arguments("contamination"),
        opt_parameters=config["genome_dereplication"]["opt_parameters"],
        min_overlap=config["genome_dereplication"]["overlap"],
    shell:
        " galah cluster "
        " --genome-fasta-directory {input.dir}"
        " --genome-fasta-extension fasta "
        " --genome-info {input.quality} "
        " --ani {params.ANI} "
        " --min-aligned-fraction {params.min_overlap} "
        " {params.opt_parameters} "
        " --min-completeness {params.min_completeness} "
        " --max-contamination {params.max_contamination} "
        " --quality-formula {params.quality_formula} "
        " --threads {threads} "
        " --output-representative-fasta-directory {output.dir} "
        " --output-cluster-definition {output.mapping_file} "
        " &> {log} "


localrules:
    rename_genomes,


# TODO: adapt rename script for galah output
checkpoint rename_genomes:
    input:
        genomes="genomes/dereplicated_genomes",
        mapping_file="genomes/clustering/allbins2genome_oldname.tsv",
    output:
        dir=directory("genomes/genomes"),
        mapfile_contigs="genomes/clustering/contig2genome.tsv",
        mapfile_genomes="genomes/clustering/old2newID.tsv",
        mapfile_bins="genomes/clustering/allbins2genome.tsv",
    params:
        rename_contigs=config["rename_mags_contigs"],
    shadow:
        "shallow"
    log:
        "logs/genomes/rename_genomes.log",
    script:
        "../scripts/rename_genomes.py"


def get_genome_dir():

    if ("genome_dir" in config) and (config["genome_dir"] is not None):

        genome_dir = config["genome_dir"]
        assert os.path.exists(genome_dir), f"{genome_dir} Doesn't exists"

        logger.info(f"Set genomes from {genome_dir}.")

        # check if genomes are present
        genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome

        if len(genomes) == 0:
            logger.error(f"No genomes found with fasta extension in {genome_dir} ")
            exit(1)

    else:
        genome_dir = "genomes/genomes"

    return genome_dir


genome_dir = get_genome_dir()


def get_all_genomes(wildcards):

    global genome_dir

    if genome_dir == "genomes/genomes":
        checkpoints.rename_genomes.get()

    # check if genomes are present
    genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome

    if len(genomes) == 0:
        logger.error(
            f"No genomes found with fasta extension in {genome_dir} "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        exit(1)

    return genomes


rule get_contig2genomes:
    input:
        genome_dir,
    output:
        "genomes/clustering/contig2genome.tsv",
    run:
        from glob import glob

        fasta_files = glob(input[0] + "/*.f*")

        with open(output[0], "w") as out_contigs:
            for fasta in fasta_files:

                bin_name, ext = os.path.splitext(os.path.split(fasta)[-1])
                # if gz remove also fasta extension
                if ext == ".gz":
                    bin_name = os.path.splitext(bin_name)[0]

                # write names of contigs in mapping file
                with open(fasta) as f:
                    for line in f:
                        if line[0] == ">":
                            header = line[1:].strip().split()[0]
                            out_contigs.write(f"{header}\t{bin_name}\n")


# alternative way to get to contigs2genomes for quantification with external genomes


ruleorder: get_contig2genomes > rename_genomes


rule run_all_checkm_lineage_wf:
    input:
        touched_output="logs/checkm_init.txt",
        dir=genome_dir,
    output:
        "genomes/checkm/completeness.tsv",
        "genomes/checkm/storage/tree/concatenated.fasta",
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        tmpdir=config["tmpdir"],
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        time=config["runtime"]["long"],
        mem=config["large_mem"],
    log:
        "logs/genomes/checkm.log",
    benchmark:
        "logs/benchmarks/checkm_lineage_wf/all_genomes.tsv"
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --tmpdir {params.tmpdir} \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.dir} \
            {params.output_dir} &> {log}
        """


### Quantification


rule build_db_genomes:
    input:
        genome_dir,
    output:
        index="ref/genome/3/summary.txt",
        fasta=temp("genomes/all_contigs.fasta"),
    threads: config.get("threads", 6)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    log:
        "logs/genomes/mapping/build_bbmap_index.log",
    shell:
        """
        cat {input}/*.fasta > {output.fasta} 2> {log}
        bbmap.sh build=3 -Xmx{resources.java_mem}G ref={output.fasta} threads={threads} local=f 2>> {log}

        """


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_MAGs:
    input:
        reads=get_quality_controlled_reads,
        ref=rules.build_db_genomes.output.index,
    output:
        sam=temp("genomes/alignments/{sample}.sam"),
        unmapped=expand(
            "genomes/alignments/unmapped/{{sample}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    params:
        input=lambda wc, input: input_params_for_bbwrap(input.reads),
        unmapped=lambda wc, output: io_params_for_tadpole(output.unmapped, "outu"),
        maxsites=config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        max_distance_between_pairs=config.get(
            "contig_max_distance_between_pairs", CONTIG_MAX_DISTANCE_BETWEEN_PAIRS
        ),
        paired_only=(
            "t"
            if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY)
            else "f"
        ),
        ambiguous="all" if CONTIG_COUNT_MULTI_MAPPED_READS else "best",
        min_id=config.get("contig_min_id", CONTIG_MIN_ID),
        maxindel=100,  # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/genomes/mapping/map_{sample}.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        """
        bbwrap.sh \
        build=3 \
        {params.input} \
        trimreaddescriptions=t \
        outm={output.sam} \
        {params.unmapped} \
        threads={threads} \
        pairlen={params.max_distance_between_pairs} \
        pairedonly={params.paired_only} \
        minid={params.min_id} \
        mdtag=t \
        xstag=fs \
        nmtag=t \
        sam=1.3 \
        ambiguous={params.ambiguous} \
        secondary=t \
        saa=f \
        maxsites={params.maxsites} \
        -Xmx{resources.java_mem}G \
        2> {log}
        """


rule pileup_MAGs:
    input:
        sam="genomes/alignments/{sample}.sam",
    output:
        # basecov=temp("genomes/alignments/{sample}_base_coverage.txt.gz"),
        # covhist=temp("genomes/alignments/{sample}_coverage_histogram.txt"),
        covstats=temp("genomes/alignments/{sample}_coverage.txt"),
        bincov=temp("genomes/alignments/{sample}_coverage_binned.txt"),
    log:
        "logs/genomes/alignments/pilup_{sample}.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        """pileup.sh in={input.sam} \
        threads={threads} \
        -Xmx{resources.java_mem}G \
        covstats={output.covstats} \
        concise=t \
        bincov={output.bincov} 2> {log}"""


rule combine_coverages_MAGs:
    input:
        covstats=expand("genomes/alignments/{sample}_coverage.txt", sample=SAMPLES),
        binned_coverage_files=expand(
            "genomes/alignments/{sample}_coverage_binned.txt", sample=SAMPLES
        ),
        contig2genome="genomes/clustering/contig2genome.tsv",
    params:
        samples=SAMPLES,
    output:
        counts="genomes/counts/raw_counts_genomes.tsv",
        binned_cov="genomes/counts/binned_coverage.tsv.gz",
        median_abund="genomes/counts/median_coverage_genomes.tsv",
    log:
        "logs/genomes/counts/combine_coverages_MAGs.log",
    threads: 1
    script:
        "../scripts/combine_coverage_MAGs.py"


# rule predict_genes_genomes:
#     input:
#         dir= genomes_dir
#     output:
#         directory("genomes/annotations/genes")
#     conda:
#         "%s/prodigal.yaml" % CONDAENV
#     log:
#         "logs/genomes/prodigal.log"
#     shadow:
#         "shallow"
#     threads:
#         config.get("threads", 1)
#     script:
#         "predict_genes_of_genomes.py"


rule predict_genes_genomes:
    input:
        os.path.join(genome_dir, "{genome}.fasta"),
    output:
        fna="genomes/annotations/genes/{genome}.fna",
        faa="genomes/annotations/genes/{genome}.faa",
        gff=temp("genomes/annotations/genes/{genome}.gff"),
    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "logs/genomes/prodigal/{genome}.txt",
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """


def get_all_genes(wildcards, extension=".faa"):
    return expand(
        "genomes/annotations/genes/{genome}{extension}",
        genome=get_all_genomes(wildcards),
        extension=extension,
    )


localrules:
    all_prodigal,


rule all_prodigal:
    input:
        get_all_genes,
    output:
        touch("genomes/annotations/genes/predicted"),
