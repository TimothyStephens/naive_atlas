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
        """
        (mkdir -p {output} && cp {input}/MAG_prokaryotic* {output}) 1>{log} 2>&1
        """


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
    conda:
        "../envs/gtdbtk.yaml"
    benchmark:
        "benchmarks/genomes/annotations/genomes/taxonomy/gtdbtk/identify.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/identify.log",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fa",
    shell:
        """
        (
        export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}"
        gtdbtk identify \\
            --genome_dir {input.genes_flag} \\
            --out_dir {params.outdir} \\
            --extension {params.extension} \\
            --tmpdir {resources.tmpdir} \\
            --cpus {threads}
        ) 1>{log[0]} 2>&1
        """


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
    conda:
        "../envs/gtdbtk.yaml"
    benchmark:
        "benchmarks/genomes/annotations/genomes/taxonomy/gtdbtk/align.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/align.log",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
    shell:
        """
        (
        export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}"
        gtdbtk align \\
            --identify_dir {params.outdir} \\
            --out_dir {params.outdir} \\
            --tmpdir {resources.tmpdir} \\
            --cpus {threads}
        ) 1>{log[0]} 2>&1
        """


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
    conda:
        "../envs/gtdbtk.yaml"
    benchmark:
        "benchmarks/genomes/annotations/genomes/taxonomy/gtdbtk/classify.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/classify.log",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fa",
    shell:
        """
        (
        export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}"
        gtdbtk classify \\
            --genome_dir {input.genome_dir} \\
            --align_dir {params.outdir} \\
            --out_dir {params.outdir} \\
            --tmpdir {resources.tmpdir} \\
            --extension {params.extension} \\
            --cpus {threads}
        ) 1>{log[0]} 2>&1
        """


localrules:
    combine_taxonomy,

rule combine_taxonomy:
    input:
        folder=f"{gtdb_dir}/classify",
    output:
        combined=f"{gtdb_dir}/gtdbtk.combined.summary.tsv",
        taxonomy="genomes/annotations/genomes/taxonomy/gtdb_taxonomy.tsv",
    log:
        "logs/genomes/annotations/genomes/taxonomy/gtdbtk/combine.log",
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
    benchmark:
        "benchmarks/genomes/annotations/genomes/tree/{msa}.tsv",
    log:
        "logs/genomes/annotations/genomes/tree/{msa}.log",
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_gtdbtk_tree", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_gtdbtk_tree", "account"),
    params:
        outdir=lambda wc, output: Path(output[0]).parent,
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        (
        export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}"
        gtdbtk infer --msa_file {input} \\
            --out_dir {params.outdir} \\
            --prefix {wildcards.msa} \\
            --cpus {threads} \\
            --tmpdir {resources.tmpdir}
        ) 1>{log} 2>&1
        """


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
    benchmark:
        "benchmarks/genomes/annotations/genomes/tree/root_tree_{msa}.tsv",
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
        database=rules.mmseqs2_download_db.output.database,
    output:
        codon=temp("genomes/annotations/{dataset}/{genome}.metaeuk.codon.fas.gz"),
        fas=temp("genomes/annotations/{dataset}/{genome}.metaeuk.fas.gz"),
        gff=temp("genomes/annotations/{dataset}/{genome}.metaeuk.gff.gz"),
        headerMap=temp("genomes/annotations/{dataset}/{genome}.metaeuk.headersMap.tsv.gz"),
        headerMap_combined=temp("genomes/annotations/{dataset}/{genome}.metaeuk_combined.headersMap.tsv.gz"),
        contig_classification=temp("genomes/annotations/{dataset}/{genome}.metaeuk_contig_classification.tsv.gz"),
        mag_classification=temp("genomes/annotations/{dataset}/{genome}.metaeuk_mag_classification.tsv.gz"),
        tax_per_contig=temp("genomes/annotations/{dataset}/{genome}.metaeuk_tax_per_contig.tsv.gz"),
        tax_per_pred=temp("genomes/annotations/{dataset}/{genome}.metaeuk_tax_per_pred.tsv.gz"),
        tax_per_contig_combined=temp("genomes/annotations/{dataset}/{genome}.metaeuk_combined_tax_per_contig.tsv.gz"),
        tax_per_pred_combined=temp("genomes/annotations/{dataset}/{genome}.metaeuk_combined_tax_per_pred.tsv.gz"),
    params:
        codon="genomes/annotations/{dataset}/{genome}.metaeuk.codon.fas",
        fas="genomes/annotations/{dataset}/{genome}.metaeuk.fas",
        gff="genomes/annotations/{dataset}/{genome}.metaeuk.gff",
        headerMap="genomes/annotations/{dataset}/{genome}.metaeuk.headersMap.tsv",
        headerMap_combined="genomes/annotations/{dataset}/{genome}.metaeuk_combined.headersMap.tsv",
        contig_classification="genomes/annotations/{dataset}/{genome}.metaeuk_contig_classification.tsv",
        mag_classification="genomes/annotations/{dataset}/{genome}.metaeuk_mag_classification.tsv",
        tax_per_contig="genomes/annotations/{dataset}/{genome}.metaeuk_tax_per_contig.tsv",
        tax_per_pred="genomes/annotations/{dataset}/{genome}.metaeuk_tax_per_pred.tsv",
        tax_per_contig_combined="genomes/annotations/{dataset}/{genome}.metaeuk_combined_tax_per_contig.tsv",
        tax_per_pred_combined="genomes/annotations/{dataset}/{genome}.metaeuk_combined_tax_per_pred.tsv",
        metaeuk_createdb=config["metaeuk_createdb"],
        metaeuk_predictexons=config["metaeuk_predictexons"],
        metaeuk_reduceredundancy=config["metaeuk_reduceredundancy"],
        metaeuk_unitesetstofasta=config["metaeuk_unitesetstofasta"],
        metaeuk_taxtocontig=config["metaeuk_taxtocontig"],
        out="genomes/annotations/{dataset}/{genome}.metaeuk",
        out_combined="genomes/annotations/{dataset}/{genome}.metaeuk_combined",
        mag_id=lambda wc: wc.genome,
        tmp="genomes/annotations/{dataset}/{genome}.tmp",
    threads: lambda wc: get_resource(wc, None, 1, "genome_annot_metaeuk", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "genome_annot_metaeuk", "account"),
    container:
        "docker://ghcr.io/soedinglab/metaeuk:7-bba0d80"
    log:
        "logs/genomes/annotations/{dataset}/{genome}.metaeuk.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/{genome}.metaeuk.tsv",
    shell:
        """
        (
        rm -fr "{params.tmp}"
        mkdir -p "{params.tmp}"
        /usr/local/bin/entrypoint createdb \\
          {input.fasta} {params.tmp}/contigDB \\
          {params.metaeuk_createdb}
        /usr/local/bin/entrypoint predictexons \\
          {params.tmp}/contigDB {input.database} {params.tmp}/callsResultDB {params.tmp}/tmp \\
          --threads {threads} {params.metaeuk_predictexons}
        /usr/local/bin/entrypoint reduceredundancy \\
          {params.tmp}/callsResultDB {params.tmp}/predsResultDB {params.tmp}/predGroupsDB \\
          --threads {threads} {params.metaeuk_reduceredundancy}
        /usr/local/bin/entrypoint unitesetstofasta \\
          {params.tmp}/contigDB {input.database} {params.tmp}/predsResultDB {params.out} \\
          --threads {threads} {params.metaeuk_unitesetstofasta}
        if [ $(grep -c '>' "{params.fas}") -gt 0 ]; then
            /usr/local/bin/entrypoint taxtocontig \\
              {params.tmp}/contigDB {params.fas} {params.headerMap} {input.database} {params.out} {params.tmp}/tmp \\
              --threads {threads} {params.metaeuk_taxtocontig}
            awk '{{OFS=FS="\\t"}}{{$1=0; print}}' {params.headerMap} > {params.headerMap_combined}
            /usr/local/bin/entrypoint taxtocontig \\
              {params.tmp}/contigDB {params.fas} {params.headerMap_combined} {input.database} {params.out_combined} {params.tmp}/tmp \\
              --threads {threads} {params.metaeuk_taxtocontig}
            # Format results
            awk 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "contig_name\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }} {{
                    print
                }}' \\
              {params.tax_per_contig} > {params.contig_classification}
            awk -v S="{params.mag_id}" 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }}NR==1{{
                    $1=S; 
                    print
                }}' \\
              {params.tax_per_contig_combined} > {params.mag_classification}
        else
            echo "WARNING: No predicted genes. Can't assess taxonomy!"
            touch "{params.headerMap_combined}"
            touch "{params.tax_per_contig}"
            touch "{params.tax_per_contig_combined}"
            # Format results
            grep '>' "{input.fasta}" | sed -e 's/>//' \\
              | awk 'BEGIN{{
                    OFS=FS="\\t"; 
                    print "contig_name\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage"
                }} {{
                    print $1"\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA"
                }}' \\
              > {params.contig_classification}
            echo -e "MAG_ID\\tmetaeuk_tax_id\\tmetaeuk_tax_rank\\tmetaeuk_name\\tmetaeuk_total_frags\\tmetaeuk_assigned_frags\\tmetaeuk_frags_agreement\\tmetaeuk_agreement_ratio\\tmetaeuk_lineage" > {params.mag_classification}
            echo -e "{params.mag_id}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" >> {params.mag_classification}
        fi
        gzip -9 {params.out}*
        rm -fr {params.tmp}
        ) 1>{log} 2>&1
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
    return(expand('genomes/annotations/{dataset}/{genome}.metaeuk_contig_classification.tsv.gz', 
                        dataset=wildcards.dataset, genome=all_genomes)
    )

def get_all_genome_metaeuk_mag_results(wildcards):
    all_genomes = get_all_genome_metaeuk(wildcards)
    return(expand('genomes/annotations/{dataset}/{genome}.metaeuk_mag_classification.tsv.gz', 
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
        "logs/genomes/annotations/{dataset}/metaeuk_combine.log",
    script:
        "../scripts/combine_metaeuk.py"





###############################
####                       ####
####        MMSEQS2        ####
####                       ####
###############################

rule genome_mmseqs2_easy_taxonomy:
    input:
        fasta="genomes/{dataset}/{genome}.fa",
        database=rules.mmseqs2_download_db.output.database,
    output:
        result_lca=temp("genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}_result_lca.tsv.gz"),
        result_report=temp("genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}_result_report.gz"),
        result_tophit_aln=temp("genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}_result_tophit_aln.gz"),
        result_tophit_report=temp("genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}_result_tophit_report.gz"),
    params:
        out="genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}_result",
        opts=config["mmseqs2_easy_taxonomy_opts"],
        mag_id=lambda wc: wc.genome,
        tmp="genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}.tmp",
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
        "logs/genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}.log",
    benchmark:
        "logs/benchmarks/genomes/annotations/{dataset}/{genome}.mmseqs2_easy_taxonomy_{database_name}.tsv",
    shell:
        """
        (
        rm -fr "{params.tmp}"
        mkdir -p "{params.tmp}"
        /usr/local/bin/entrypoint easy-taxonomy \\
          {input.fasta} {input.database} \\
          {params.out} {params.tmp} \\
          {params.opts} \\
          --threads {threads} \\
          --split-memory-limit {resources.mem_gb}G \\
        && gzip -9 {params.out}* \\
        && rm -fr {params.tmp}
        ) 1>{log} 2>&1
        """


def get_all_genome_mmseqs2_easy_taxonomy_results(wildcards):
    if wildcards.dataset == "genomes":
        all_genomes = get_all_genomes(wildcards)
    else:
        all_genomes = get_all_unbinned(wildcards)
    
    return {
        "result_lca": expand(
            rules.genome_mmseqs2_easy_taxonomy.output.result_lca,
            dataset=wildcards.dataset, genome=all_genomes, database_name=wildcards.database_name
        ),
        "result_tophit_aln": expand(
            rules.genome_mmseqs2_easy_taxonomy.output.result_tophit_aln,
            dataset=wildcards.dataset, genome=all_genomes, database_name=wildcards.database_name
        ),
        "result_tophit_report": expand(
            rules.genome_mmseqs2_easy_taxonomy.output.result_tophit_report,
            dataset=wildcards.dataset, genome=all_genomes, database_name=wildcards.database_name
        )
    }

localrules:
    all_genome_mmseqs2_easy_taxonomy,

rule all_genome_mmseqs2_easy_taxonomy:
    input:
        unpack(get_all_genome_mmseqs2_easy_taxonomy_results),
    output:
        result_lca = "genomes/annotations/{dataset}/mmseqs2_easy_taxonomy_{database_name}_result_lca.tsv.gz",
        result_tophit_aln = "genomes/annotations/{dataset}/mmseqs2_easy_taxonomy_{database_name}_result_tophit_aln.tsv.gz",
        result_tophit_report = "genomes/annotations/{dataset}/mmseqs2_easy_taxonomy_{database_name}_result_tophit_report.tsv.gz"
    log:
        "logs/genomes/annotations/{dataset}/mmseqs2_easy_taxonomy_combine_{database_name}.log"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account")
    run:
        try:
            import traceback
            import pandas as pd
            
            def combine_and_save(file_list, out_file, col_names):
                tables = [
                    pd.read_csv(f, sep="\t", header=None, names=col_names)
                    for f in file_list
                ]
                combined = pd.concat(tables, axis=0, ignore_index=True)
                combined.to_csv(out_file, sep='\t', index=False, header=True)
            
            combine_and_save(input.result_lca, output.result_lca,
                ["query", "taxid", "rank", "name", "retained_taxa", "agreement", "evalue"]
            )
            combine_and_save(input.result_tophit_aln, output.result_tophit_aln,
                ["query", "target", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]
            )
            combine_and_save(input.result_tophit_report, output.result_tophit_report,
                ["percent_reads", "clade_reads", "taxon_reads", "rank", "taxid", "name"]
            )
        
        except Exception as e:
            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)
            raise e

