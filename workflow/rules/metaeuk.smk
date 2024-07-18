

localrules:
    MetaEuk_download,


rule metaeuk_download:
    output:
        dbdir=directory(f"{DBDIR}/MetaEuk"),
        database=os.path.join(f"{DBDIR}/MetaEuk", config["metaeuk_database"]),
    params:
        metaeuk_database=config["metaeuk_database"],
    threads: config["large_threads"]
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    log:
        "logs/MetaEuk/download_MetaEuk_database.log",
    benchmark:
        "logs/benchmarks/MetaEuk/download_MetaEuk_database.tsv"
    conda:
        "../envs/metaeuk.yaml"
    shell:
        "metaeuk databases {params.metaeuk_database} {output.database} {output.dbdir}/tmp "
        " --compressed 1 "
        " --threads {threads} "
        " &> {log}"


rule metaeuk_annotation:
    input:
        fasta="genomes/genomes/{genome}.fasta",
        database=rules.metaeuk_download.output.database,
    output:
        codon="genomes/annotations/metaeuk/{genome}.fasta.metaeuk.codon.fas",
        fas="genomes/annotations/metaeuk/{genome}.fasta.metaeuk.fas",
        gff="genomes/annotations/metaeuk/{genome}.fasta.metaeuk.gff",
        headerMap="genomes/annotations/metaeuk/{genome}.fasta.metaeuk.headersMap.tsv",
        tax_per_contig="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_tax_per_contig.tsv",
        tax_per_pred="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_tax_per_pred.tsv",
        tmp=temp(directory("genomes/annotations/metaeuk/{genome}.fasta.metaeuk.tmp")),
    params:
        metaeuk_createdb=config["metaeuk_createdb"],
        metaeuk_predictexons=config["metaeuk_predictexons"],
        metaeuk_reduceredundancy=config["metaeuk_reduceredundancy"],
        metaeuk_unitesetstofasta=config["metaeuk_unitesetstofasta"],
        metaeuk_taxtocontig=config["metaeuk_taxtocontig"],
        out="genomes/annotations/metaeuk/{genome}.fasta.metaeuk"
    threads: config["large_threads"]
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/metaeuk.yaml"
    log:
        "logs/MetaEuk/run_MetaEuk/{genome}.log",
    benchmark:
        "logs/benchmarks/MetaEuk/run_MetaEuk/{genome}.tsv"
    shell:
        "("
        "mkdir -p {output.tmp}; "
        "metaeuk createdb"
        "  {input.fasta} {output.tmp}/contigDB"
        "  {params.metaeuk_createdb}; "
        "metaeuk predictexons"
        "  {output.tmp}/contigDB {input.database} {output.tmp}/callsResultDB {output.tmp}/tmp"
        "  --threads {threads} {params.metaeuk_predictexons}; "
        "metaeuk reduceredundancy"
        "  {output.tmp}/callsResultDB {output.tmp}/predsResultDB {output.tmp}/predGroupsDB"
        "  --threads {threads} {params.metaeuk_reduceredundancy}; "
        "metaeuk unitesetstofasta"
        "  {output.tmp}/contigDB {input.database} {output.tmp}/predsResultDB {params.out}"
        "  --threads {threads} {params.metaeuk_unitesetstofasta}; "
        "metaeuk taxtocontig"
        "  {output.tmp}/contigDB {output.fas} {output.headerMap} {input.database} {params.out} {output.tmp}/tmp"
        "  --threads {threads} {params.metaeuk_taxtocontig}; "
        ") &> {log}"


def get_all_metaeuk(wildcards):
    all_genomes = get_all_genomes(wildcards)
    all_unbinned_genomes = get_all_unbinned(wildcards)
    all_genomes.extend(['unbinned/'+l for l in all_unbinned_genomes])
    
    return (
        expand(rules.metaeuk_annotation.output.codon, genome=all_genomes) +
        expand(rules.metaeuk_annotation.output.fas, genome=all_genomes) +
        expand(rules.metaeuk_annotation.output.gff, genome=all_genomes) +
        expand(rules.metaeuk_annotation.output.headerMap, genome=all_genomes) +
        expand(rules.metaeuk_annotation.output.tax_per_contig, genome=all_genomes) +
        expand(rules.metaeuk_annotation.output.tax_per_pred, genome=all_genomes)
    )

localrules:
    all_metaeuk,


rule all_metaeuk:
    input:
        get_all_metaeuk,
    output:
        touch("genomes/annotations/metaeuk/finished"),


