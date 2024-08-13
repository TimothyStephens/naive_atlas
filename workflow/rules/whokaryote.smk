


rule whokaryote_annotation:
    input:
        fasta="genomes/genomes/{genome}.fasta",
    output:
        pred="genomes/annotations/whokaryote/{genome}.fasta.whokaryote_predictions.tsv",
        predT="genomes/annotations/whokaryote/{genome}.fasta.whokaryote/whokaryote_predictions_T.tsv",
    params:
        minsize=config["whokaryote_min_contig_length"],
        out="genomes/annotations/whokaryote/{genome}.fasta.whokaryote",
        mag_id=lambda wildcards: wildcards.genome,
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/whokaryote.yaml"
    log:
        "logs/genomes/annotations/whokaryote/{genome}.log",
    benchmark:
        "logs/benchmarks/whokaryote/{genome}.tsv"
    shell:
        """
        (
        whokaryote.py \
          --contigs {input.fasta} \
          --outdir {params.out} \
          --minsize {params.minsize} \
          --threads {threads}
        echo -e "contig_id\\twhokaryote_predicted" > {output.pred}
        awk 'NR>1' {output.predT} >> {output.pred}
        ) &> {log}
        """



def get_all_whokaryote(wildcards):
    all_genomes = get_all_genomes(wildcards)
    all_unbinned_genomes = get_all_unbinned(wildcards)
    all_genomes.extend(['unbinned/'+l for l in all_unbinned_genomes])
    return all_genomes

def get_all_whokaryote_contigs(wildcards):
    all_genomes = get_all_whokaryote(wildcards)
    return( expand('genomes/genomes/{genome}.fasta', genome=all_genomes) )

def get_all_whokaryote_results(wildcards):
    all_genomes = get_all_whokaryote(wildcards)
    return( expand('genomes/annotations/whokaryote/{genome}.fasta.whokaryote_predictions.tsv', genome=all_genomes) )

rule combine_whokaryote:
    input:
        contig_fasta_files=get_all_whokaryote_contigs,
        results_files=get_all_whokaryote_results,
    output:
        output_table="genomes/annotations/whokaryote_predictions.tsv",
    params:
        genomes=get_all_whokaryote,
    log:
        "logs/genomes/annotations/whokaryote/combine.log",
    script:
        "../scripts/combine_whokaryote.py"



localrules:
    all_whokaryote,


rule all_whokaryote:
    input:
        rules.combine_whokaryote.output.output_table,
    output:
        touch("genomes/annotations/whokaryote/finished"),


