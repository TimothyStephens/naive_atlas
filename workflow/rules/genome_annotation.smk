from pathlib import Path




###############################
####                       ####
####         GTDBTK        ####
####                       ####
###############################

gtdb_dir = "genomes/annotations/genomes/taxonomy/gtdb"

localrules:
    copy_prokaryotic_genomes,

rule copy_prokaryotic_genomes:
    input:
        "genomes/genomes",
    output:
        directory("tmp/gtdbtk"),
    log:
        "logs/genomes/annotations/genomes/copy_prokaryotic_genomes.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        "mkdir -p {output} && cp {input}/MAG_prokaryotic* {output}"


rule identify:
    input:
        flag=rules.gtdb_extract.output,
        #flag=rules.extract_gtdb.output,
        genes_flag=rules.copy_prokaryotic_genomes.output,
    output:
        directory(f"{gtdb_dir}/identify"),
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "account"),
        tmpdir=config.get("tmpdir", "/tmp"),
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
        "--tmpdir {resources.tmpdir} "
        "--cpus {threads} &> {log[0]}"


checkpoint align:
    input:
        f"{gtdb_dir}/identify",
    output:
        directory(f"{gtdb_dir}/align"),
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "account"),
        tmpdir=config.get("tmpdir", "/tmp"),
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
        "--tmpdir {resources.tmpdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=rules.copy_prokaryotic_genomes.output,
    output:
        directory(f"{gtdb_dir}/classify"),
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk", "account"),
        tmpdir=config.get("tmpdir", "/tmp"),
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fa",
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--tmpdir {resources.tmpdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


localrules:
    combine_taxonomy,

rule combine_taxonomy:
    input:
        folder=f"{gtdb_dir}/classify",
    output:
        combined=f"{gtdb_dir}/gtdbtk.combined.summary.tsv",
        taxonomy="genomes/annotations/genomes/taxonomy/gtdb_taxonomy.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/combine.txt",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
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
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk_tree", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "account"),
        tmpdir=config.get("tmpdir", "/tmp"),
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
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk_tree", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "account"),
    log:
        "logs/genomes/annotations/genomes/tree/root_tree_{msa}.log",
    script:
        "../scripts/root_tree.py"


def all_gtdb_trees_input(wildcards):
    dir = checkpoints.align.get().output[0]

    domains = glob_wildcards(f"{dir}/gtdbtk.{{domain}}.user_msa.fasta.gz").domain

    return expand("genomes/annotations/genomes/tree/gtdbtk.{domain}.nwk", domain=domains)


localrules:
    all_gtdb_trees,

rule all_gtdb_trees:
    input:
        all_gtdb_trees_input,
    output:
        touch("genomes/annotations/genomes/tree/finished_gtdb_trees"),
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),





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
        codon="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.codon.fas.gz",
        fas="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.fas.gz",
        gff="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.gff.gz",
        headerMap="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.headersMap.tsv.gz",
        headerMap_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined.headersMap.tsv.gz",
        contig_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_contig_classification.tsv.gz",
        mag_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_mag_classification.tsv.gz",
        tmp=temp(directory("genomes/annotations/{dataset}/metaeuk/{genome}.fa.tmp")),
    params:
        codon="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.codon.fas",
        fas="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.fas",
        gff="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.gff",
        headerMap="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk.headersMap.tsv",
        headerMap_combined="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_combined.headersMap.tsv",
        contig_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_contig_classification.tsv",
        mag_classification="genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_mag_classification.tsv",
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
        mag_id=lambda wc: wc.genome,
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_metaeuk", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "account"),
    container:
        "docker://ghcr.io/soedinglab/metaeuk:7-bba0d80"
    log:
        "logs/genomes/annotations/{dataset}/metaeuk/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/metaeuk/{genome}.tsv"
    shell:
        """
        (
        mkdir -p {output.tmp}
        /usr/local/bin/entrypoint createdb \
          {input.fasta} {output.tmp}/contigDB \
          {params.metaeuk_createdb}
        /usr/local/bin/entrypoint predictexons \
          {output.tmp}/contigDB {input.database} {output.tmp}/callsResultDB {output.tmp}/tmp \
          --threads {threads} {params.metaeuk_predictexons}
        /usr/local/bin/entrypoint reduceredundancy \
          {output.tmp}/callsResultDB {output.tmp}/predsResultDB {output.tmp}/predGroupsDB \
          --threads {threads} {params.metaeuk_reduceredundancy}
        /usr/local/bin/entrypoint unitesetstofasta \
          {output.tmp}/contigDB {input.database} {output.tmp}/predsResultDB {params.out} \
          --threads {threads} {params.metaeuk_unitesetstofasta}
        if [ $(grep -c '>' "{params.fas}") -gt 0 ]; then
            /usr/local/bin/entrypoint taxtocontig \
              {output.tmp}/contigDB {params.fas} {params.headerMap} {input.database} {params.out} {output.tmp}/tmp \
              --threads {threads} {params.metaeuk_taxtocontig}
            awk '{{OFS=FS="\\t"}}{{$1=0; print}}' {params.headerMap} > {params.headerMap_combined}
            /usr/local/bin/entrypoint taxtocontig \
              {output.tmp}/contigDB {params.fas} {params.headerMap_combined} {input.database} {params.out_combined} {output.tmp}/tmp \
              --threads {threads} {params.metaeuk_taxtocontig}
            # Format results
            awk 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "contig_name\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }} {{
                    print
                }}' \
              {params.tax_per_contig} > {params.contig_classification}
            awk -v S="{params.mag_id}" 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }}NR==1{{
                    $1=S; 
                    print
                }}' \
              {params.tax_per_contig_combined} > {params.mag_classification}
        else
            echo "WARNING: No predicted genes. Can't assess taxonomy!"
            touch "{params.headerMap_combined}"
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
              > {params.contig_classification}
            echo -e "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage" > {params.mag_classification}
            echo -e "{params.mag_id}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" >> {params.mag_classification}
        fi
        gzip -9 {params.out}*
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
    return(expand('genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_contig_classification.tsv.gz', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

def get_all_genome_metaeuk_mag_results(wildcards):
    all_genomes = get_all_genome_metaeuk(wildcards)
    return(expand('genomes/annotations/{dataset}/metaeuk/{genome}.fa.metaeuk_mag_classification.tsv.gz', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

rule combine_genome_metaeuk:
    input:
        contig_fasta_files=get_all_genome_metaeuk_contigs,
        contig_results_files=get_all_genome_metaeuk_contig_results,
        mag_results_files=get_all_genome_metaeuk_mag_results,
    output:
        contig_output_table="genomes/annotations/{dataset}/metaeuk_contig_predictions.tsv.gz",
        mag_output_table="genomes/annotations/{dataset}/metaeuk_mag_predictions.tsv.gz",
    params:
        genomes=get_all_genome_metaeuk,
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
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
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),





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
        result_lca="genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy_result_lca.tsv.gz",
        result_report="genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy_result_report.gz",
        result_tophit_aln="genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy_result_tophit_aln.gz",
        result_tophit_report="genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy_result_tophit_report.gz",
        tmp=temp(directory("genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy.tmp")),
    params:
        out="genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.easy_taxonomy_result",
        mmseqs2_easy_taxonomy=config["mmseqs2_easy_taxonomy"],
        mag_id=lambda wc: wc.genome,
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_mmseqs2_easy_taxonomy", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_mmseqs2_easy_taxonomy", "mem_mb"),
        mem_gb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_mmseqs2_easy_taxonomy", "mem_gb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_mmseqs2_easy_taxonomy", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_mmseqs2_easy_taxonomy", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_mmseqs2_easy_taxonomy", "account"),
    container:
        "docker://ghcr.io/soedinglab/mmseqs2:18-8cc5c"
    log:
        "logs/genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/{genome}.tsv"
    shell:
        """
        (
        mkdir -p {output.tmp}
        /usr/local/bin/entrypoint easy-taxonomy \
          {input.fasta} {input.database} \
          {params.out} {output.tmp} \
          {params.mmseqs2_easy_taxonomy} \
          --threads {threads} \
          --split-memory-limit {resources.mem_gb}G \
        && gzip -9 {params.out}*	
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
        touch("genomes/annotations/{dataset}/mmseqs2_easy_taxonomy/finished"),
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),


