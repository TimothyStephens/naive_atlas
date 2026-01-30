




###############################
####                       ####
####         GTDBTK        ####
####                       ####
###############################

gtdb_dir = "genomes/annotations/genomes/taxonomy/gtdb"

rule copy_prokaryotic_genomes:
    input:
        "genomes/genomes",
    output:
        directory("tmp/gtdbtk"),
    log:
        "logs/genomes/annotations/genomes/copy_prokaryotic_genomes.log",
    shell:
        "mkdir -p {output} && cp {input}/MAG_prokaryotic* {output}"


rule identify:
    input:
        flag=rules.gtdb_extract.output,
        #flag=rules.extract_gtdb.output,
        genes_flag=rules.copy_prokaryotic_genomes.output,
    output:
        directory(f"{gtdb_dir}/identify"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["large_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/identify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fa",
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk identify "
        "--genome_dir {input.genes_flag} "
        "--out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


checkpoint align:
    input:
        f"{gtdb_dir}/identify",
    output:
        directory(f"{gtdb_dir}/align"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["large_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/align.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=rules.copy_prokaryotic_genomes.output,
    output:
        directory(f"{gtdb_dir}/classify"),
    threads: config["simplejob_threads"]  #pplacer needs much memory for not many threads
    resources:
        mem=config["large_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fa",
        mashdir=Path(GTDBTK_DATA_PATH) / "mash_db",
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--tmpdir {resources.tmpdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


rule combine_taxonomy:
    input:
        folder=f"{gtdb_dir}/classify",
    output:
        combined=f"{gtdb_dir}/gtdbtk.combined.summary.tsv",
        taxonomy="genomes/annotations/genomes/taxonomy/gtdb_taxonomy.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/combine.txt",
    script:
        "../scripts/combine_taxonomy.py"


rule build_tree:
    input:
        f"{gtdb_dir}/align/{{msa}}.user_msa.fasta.gz",
    output:
        temp("genomes/annotations/genomes/taxonomy/gtdb/{msa}.unrooted.tree"),
    log:
        "logs/genomes/annotations/genomes/tree/{msa}.log",
        "logs/genomes/annotations/genomes/tree/{msa}.err",
    threads: max(config["simplejob_threads"], 3)
    params:
        outdir=lambda wc, output: Path(output[0]).parent,
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk infer --msa_file {input} "
        " --out_dir {params.outdir} "
        " --prefix {wildcards.msa} "
        " --cpus {threads} "
        "--tmpdir {resources.tmpdir} > {log[0]} 2> {log[1]}"


localrules:
    root_tree,


rule root_tree:
    input:
        tree=rules.build_tree.output[0],
    wildcard_constraints:
        msa="((?!unrooted).)*",
    output:
        tree="genomes/annotations/genomes/tree/{msa}.nwk",
    conda:
        "../envs/tree.yaml"
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    log:
        "logs/genomes/annotations/genomes/tree/root_tree_{msa}.log",
    script:
        "../scripts/root_tree.py"


def all_gtdb_trees_input(wildcards):
    dir = checkpoints.align.get().output[0]

    domains = glob_wildcards(f"{dir}/gtdbtk.{{domain}}.user_msa.fasta.gz").domain

    return expand("genomes/annotations/genomes/tree/gtdbtk.{domain}.nwk", domain=domains)


rule all_gtdb_trees:
    input:
        all_gtdb_trees_input,
    output:
        touch("genomes/annotations/genomes/tree/finished_gtdb_trees"),





###############################
####                       ####
####         DRAM          ####
####                       ####
###############################

DBDIR = config["database_dir"]

def get_dram_config(wildcards):
    old_dram_path = f"{DBDIR}/Dram"
    if Path(old_dram_path).exists():
        logger.error(
            f"Detected an old database for DRAM in {old_dram_path}. You can delete it."
        )

    return config.get("dram_config_file", f"{DBDIR}/DRAM/DRAM.config")


rule genome_DRAM_annotate:
    input:
        fasta="genomes/{dataset}/{genome}.fa",
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/{dataset}/dram/intermediate_files/{genome}"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    params:
        extra=config.get("dram_extra", ""),
        min_contig_size=config.get("minimum_contig_length", "1000"),
    log:
        "logs/annotations/{dataset}/dram/run_dram/{genome}.log",
    benchmark:
        "logs/benchmarks/annotations/{dataset}/dram/run_dram/{genome}.tsv"
    shell:
        " DRAM.py annotate "
        " --config_loc {input.config} "
        " --input_fasta {input.fasta}"
        " --output_dir {output.outdir} "
        " --threads {threads} "
        " --min_contig_size {params.min_contig_size} "
        " {params.extra} "
        " --verbose &> {log}"
        #" --gtdb_taxonomy {input.gtdb_dir}/{params.gtdb_file} "
        #" --checkm_quality {input.checkm} "


def get_all_genome_dram(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return expand(rules.genome_DRAM_annotate.output.outdir,
            dataset=wildcards.dataset, genome=all_genomes)


localrules:
    concat_annotations,

rule concat_annotations:
    input:
        get_all_genome_dram,
    output:
        "genomes/annotations/{dataset}/dram/annotations.tsv",
    resources:
        time=config["simplejob_runtime"],
    run:
        from utils import io

        for i, annotation_file in enumerate(["annotations.tsv"]):
            input_files = [
                os.path.join(dram_folder, annotation_file) for dram_folder in input
            ]

            io.pandas_concat(
                input_files, output[i], sep="\t", index_col=0, axis=0, disk_based=True
            )


rule genome_DRAM_destill:
    input:
        rules.concat_annotations.output,
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/{dataset}/dram/distil"),
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/annotations/{dataset}/dram/distil.log",
    shell:
        " DRAM.py distill "
        " --config_loc {input.config} "
        " --input_file {input[0]}"
        " --output_dir {output} "
        "  &> {log}"


rule get_all_genome_modules:
    input:
        annotations="genomes/annotations/{dataset}/dram/annotations.tsv",
        config=get_dram_config,
    output:
        "genomes/annotations/{dataset}/dram/kegg_modules.tsv",
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/annotations/{dataset}/dram/get_all_modules.log",
    script:
        "../scripts/DRAM_get_all_modules.py"


rule dram:
    input:
        "genomes/annotations/{dataset}/dram/distil",
        "genomes/annotations/{dataset}/dram/kegg_modules.tsv",
    output:
        touch("genomes/annotations/{dataset}/dram/finished"),





###############################
####                       ####
####        MetaEuk        ####
####                       ####
###############################

rule genome_metaeuk_annotation:
    input:
        fasta="genomes/{dataset}/{genome}.fa",
        database=rules.mmseqs2_download.output.database,
    output:
        codon="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.codon.fas",
        fas="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.fas",
        gff="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.gff",
        headerMap="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.headersMap.tsv",
        headerMap_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined.headersMap.tsv",
        contig_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_contig_classification.tsv",
        mag_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_mag_classification.tsv",
        tmp=temp(directory("genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.tmp")),
    params:
        tax_per_contig="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_tax_per_contig.tsv",
        tax_per_pred="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_tax_per_pred.tsv",
        tax_per_contig_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined_tax_per_contig.tsv",
        tax_per_pred_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined_tax_per_pred.tsv",
        metaeuk_createdb=config["metaeuk_createdb"],
        metaeuk_predictexons=config["metaeuk_predictexons"],
        metaeuk_reduceredundancy=config["metaeuk_reduceredundancy"],
        metaeuk_unitesetstofasta=config["metaeuk_unitesetstofasta"],
        metaeuk_taxtocontig=config["metaeuk_taxtocontig"],
        out="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk",
        out_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined",
        mag_id=lambda wildcards: wildcards.genome,
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/metaeuk.yaml"
    log:
        "logs/genomes/annotations/{dataset}/metaeuk/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/metaeuk/{genome}.tsv"
    shell:
        """
        (
        mkdir -p {output.tmp}
        metaeuk createdb \
          {input.fasta} {output.tmp}/contigDB \
          {params.metaeuk_createdb}
        metaeuk predictexons \
          {output.tmp}/contigDB {input.database} {output.tmp}/callsResultDB {output.tmp}/tmp \
          --threads {threads} {params.metaeuk_predictexons}
        metaeuk reduceredundancy \
          {output.tmp}/callsResultDB {output.tmp}/predsResultDB {output.tmp}/predGroupsDB \
          --threads {threads} {params.metaeuk_reduceredundancy}
        metaeuk unitesetstofasta \
          {output.tmp}/contigDB {input.database} {output.tmp}/predsResultDB {params.out} \
          --threads {threads} {params.metaeuk_unitesetstofasta}
        if [ $(grep -c '>' "{output.fas}") -gt 0 ]; then
            metaeuk taxtocontig \
              {output.tmp}/contigDB {output.fas} {output.headerMap} {input.database} {params.out} {output.tmp}/tmp \
              --threads {threads} {params.metaeuk_taxtocontig}
            awk '{{OFS=FS="\\t"}}{{$1=0; print}}' {output.headerMap} > {output.headerMap_combined}
            metaeuk taxtocontig \
              {output.tmp}/contigDB {output.fas} {output.headerMap_combined} {input.database} {params.out_combined} {output.tmp}/tmp \
              --threads {threads} {params.metaeuk_taxtocontig}
            # Format results
            awk 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "contig_name\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }} {{
                    print
                }}' \
              {params.tax_per_contig} > {output.contig_classification}
            awk -v S="{params.mag_id}" 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }}NR==1{{
                    $1=S; 
                    print
                }}' \
              {params.tax_per_contig_combined} > {output.mag_classification}
        else
            echo "WARNING: No predicted genes. Can't assess taxonomy!"
            touch "{output.headerMap_combined}"
            touch "{params.tax_per_contig}"
            touch "{params.tax_per_contig_combined}"
            # Format results
            grep '>' "{input.fasta}" | sed -e 's/>//' \
              | awk 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "contig_name\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }} {{
                    print $1"\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA"
                }}' \
              > {output.contig_classification}
            echo -e "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage" > {output.mag_classification}
            echo -e "{params.mag_id}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" >> {output.mag_classification}
        fi
        ) &> {log}
        """



def get_all_genome_metaeuk(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    return all_genomes

def get_all_genome_metaeuk_contigs(wildcards):
    all_genomes = get_all_genome_metaeuk(wildcards)
    return(expand('genomes/{dataset}/{genome}.fa', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

def get_all_genome_metaeuk_contig_results(wildcards):
    all_genomes = get_all_genome_metaeuk(wildcards)
    return(expand('genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_contig_classification.tsv', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

def get_all_genome_metaeuk_mag_results(wildcards):
    all_genomes = get_all_genome_metaeuk(wildcards)
    return(expand('genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_mag_classification.tsv', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

rule combine_genome_metaeuk:
    input:
        contig_fasta_files=get_all_genome_metaeuk_contigs,
        contig_results_files=get_all_genome_metaeuk_contig_results,
        mag_results_files=get_all_genome_metaeuk_mag_results,
    output:
        contig_output_table="genomes/annotations/{dataset}/metaeuk_contig_predictions.tsv",
        mag_output_table="genomes/annotations/{dataset}/metaeuk_mag_predictions.tsv",
    params:
        genomes=get_all_genome_metaeuk,
    log:
        "logs/genomes/annotations/{dataset}/metaeuk/combine.log",
    script:
        "../scripts/combine_metaeuk.py"



localrules:
    all_genome_metaeuk,

rule all_genome_metaeuk:
    input:
        rules.combine_genome_metaeuk.output.contig_output_table,
        rules.combine_genome_metaeuk.output.mag_output_table,
    output:
        touch("genomes/annotations/{dataset}/metaeuk/finished"),





###############################
####                       ####
####        MMSEQS2        ####
####                       ####
###############################

rule genome_mmseqs2_easy_taxonomy:
    input:
        fasta="genomes/{dataset}/{genome}.fa",
        database=rules.mmseqs2_download.output.database,
    output:
        result_lca="genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy_result_lca.tsv",
        result_report="genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy_result_report",
        result_tophit_aln="genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy_result_tophit_aln",
        result_tophit_report="genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy_result_tophit_report",
        tmp=temp(directory("genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy.tmp")),
    params:
        out="genomes/annotations/{dataset}/mmseqs2/{genome}.easy_taxonomy_result",
        mmseqs2_easy_taxonomy=config["mmseqs2_easy_taxonomy"],
        mag_id=lambda wildcards: wildcards.genome,
        mem=int(config["mmseqs2_memory"]*0.8),
    threads: config["simplejob_threads"]
    resources:
        mem=config["mmseqs2_memory"],
        time=config["simplejob_runtime"],
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "logs/genomes/annotations/{dataset}/mmseqs2/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/mmseqs2/{genome}.tsv"
    shell:
        """
        (
        mkdir -p {output.tmp}
        mmseqs easy-taxonomy \
          {input.fasta} {input.database} \
          {params.out} {output.tmp} \
          {params.mmseqs2_easy_taxonomy} \
          --threads {threads} \
          --split-memory-limit {params.mem}G
        ) &> {log}
        """



def get_all_genome_mmseqs2_easy_taxonomy_results(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    
    return(expand(rules.genome_mmseqs2_easy_taxonomy.output.result_report, 
                    dataset=wildcards.dataset, genome=all_genomes)
    )



localrules:
    all_genome_mmseqs2_easy_taxonomy,

rule all_genome_mmseqs2_easy_taxonomy:
    input:
        get_all_genome_mmseqs2_easy_taxonomy_results,
    output:
        touch("genomes/annotations/{dataset}/mmseqs2/easy_taxonomy_finished"),




