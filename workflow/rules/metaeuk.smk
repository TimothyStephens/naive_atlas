

localrules:
    metaeuk_download,


rule metaeuk_download:
    output:
        dbdir=directory(f"{DBDIR}/MetaEuk"),
        database=os.path.join(f"{DBDIR}/MetaEuk", config["metaeuk_database_name"]),
    params:
        metaeuk_database=config["metaeuk_database"],
    threads: config["large_threads"]
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    log:
        "logs/genomes/annotations/metaeuk/download_MetaEuk_database.log",
    benchmark:
        "logs/benchmarks/metaeuk/download_MetaEuk_database.tsv"
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
        headerMap_combined="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_combined.headersMap.tsv",
        contig_classification="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_contig_classification.tsv",
        mag_classification="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_mag_classification.tsv",
        tmp=temp(directory("genomes/annotations/metaeuk/{genome}.fasta.metaeuk.tmp")),
    params:
        tax_per_contig="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_tax_per_contig.tsv",
        tax_per_pred="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_tax_per_pred.tsv",
        tax_per_contig_combined="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_combined_tax_per_contig.tsv",
        tax_per_pred_combined="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_combined_tax_per_pred.tsv",
        metaeuk_createdb=config["metaeuk_createdb"],
        metaeuk_predictexons=config["metaeuk_predictexons"],
        metaeuk_reduceredundancy=config["metaeuk_reduceredundancy"],
        metaeuk_unitesetstofasta=config["metaeuk_unitesetstofasta"],
        metaeuk_taxtocontig=config["metaeuk_taxtocontig"],
        out="genomes/annotations/metaeuk/{genome}.fasta.metaeuk",
        out_combined="genomes/annotations/metaeuk/{genome}.fasta.metaeuk_combined",
        mag_id=lambda wildcards: wildcards.genome,
    threads: config["large_threads"]
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/metaeuk.yaml"
    log:
        "logs/genomes/annotations/metaeuk/{genome}.log",
    benchmark:
        "logs/benchmarks/metaeuk/{genome}.tsv"
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
        if [ -s {output.headerMap} ];
        then
            metaeuk unitesetstofasta \
              {output.tmp}/contigDB {input.database} {output.tmp}/predsResultDB {params.out} \
              --threads {threads} {params.metaeuk_unitesetstofasta}
        else
            touch 
        fi
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



def get_all_metaeuk(wildcards):
    all_genomes = get_all_genomes(wildcards)
    all_unbinned_genomes = get_all_unbinned(wildcards)
    all_genomes.extend(['unbinned/'+l for l in all_unbinned_genomes])
    return all_genomes

def get_all_metaeuk_contigs(wildcards):
    all_genomes = get_all_metaeuk(wildcards)
    return( expand('genomes/genomes/{genome}.fasta', genome=all_genomes) )

def get_all_metaeuk_contig_results(wildcards):
    all_genomes = get_all_metaeuk(wildcards)
    return( expand('genomes/annotations/metaeuk/{genome}.fasta.metaeuk_contig_classification.tsv', genome=all_genomes) )

def get_all_metaeuk_mag_results(wildcards):
    all_genomes = get_all_metaeuk(wildcards)
    return( expand('genomes/annotations/metaeuk/{genome}.fasta.metaeuk_mag_classification.tsv', genome=all_genomes) )


rule combine_metaeuk:
    input:
        contig_fasta_files=get_all_metaeuk_contigs,
        contig_results_files=get_all_metaeuk_contig_results,
        mag_results_files=get_all_metaeuk_mag_results,
    output:
        contig_output_table="genomes/annotations/metaeuk_contig_predictions.tsv",
        mag_output_table="genomes/annotations/metaeuk_mag_predictions.tsv",
    params:
        genomes=get_all_metaeuk,
    log:
        "logs/genomes/annotations/metaeuk/combine.log",
    script:
        "../scripts/combine_metaeuk.py"



localrules:
    all_metaeuk,


rule all_metaeuk:
    input:
        rules.combine_metaeuk.output.contig_output_table,
        rules.combine_metaeuk.output.mag_output_table,
    output:
        touch("genomes/annotations/metaeuk/finished"),


