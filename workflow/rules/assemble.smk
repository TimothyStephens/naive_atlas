import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings
from copy import deepcopy


def get_preprocessing_steps(config):
    preprocessing_steps = ["QC"]
    if config.get("normalize_reads_before_assembly", False):
        preprocessing_steps.append("normalized")

    if config.get("error_correction_before_assembly", True):
        preprocessing_steps.append("errorcorr")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)

#
rule normalize_reads:
    input:
        reads=expand(
            "{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    output:
        reads=temp(
            expand(
                "{{sample}}/assembly/reads/{{previous_steps}}.normalized_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
        histin="{sample}/assembly/normalization/histogram_{previous_steps}_before_normalization.tsv.gz",
        histout=(
            "{sample}/assembly/normalization/histogram_{previous_steps}_after.tsv.gz"
        ),
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input.reads),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, key="out"),
        pairs=lambda wc, input, output: " ".join(",".join(x) for x in list(zip(input.reads, output.reads))),
        run_step=lambda wc: "t" if check_bool(wc, "Normalize_reads") else "f",
        k=config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
        target=config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
        mindepth=config["normalization_minimum_kmer_depth"],
    log:
        "{sample}/logs/assembly/pre_process/normalization_{previous_steps}.log",
    benchmark:
        "logs/benchmarks/assembly/pre_process/normalization/{sample}_{previous_steps}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ( 
        if [ '{params.run_step}' == 't' ]; then
            bbnorm.sh {params.inputs} \
                {params.outputs} \
                tmpdir={resources.tmpdir} \
                tossbadreads=t \
                hist={output.histin} \
                histout={output.histout} \
                mindepth={params.mindepth} \
                k={params.k} \
                target={params.target} \
                prefilter=t \
                threads={threads} \
                -Xmx{resources.java_mem}G
        else
            echo 'NOTE: Skipping step and hard linking files instead.'
            mkdir -p "{output.reads}"
            touch "{output.stats}"
            for group in {params.pairs};
            do
                IFS=',' read -ra items <<< "$group"
                ln ${{items[0]}} ${{items[1]}}
            done
        fi
        ) 1>{log} 2>&1
        """


rule error_correction:
    input:
        reads=expand(
            "{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    output:
        reads=temp(
            expand(
                "{{sample}}/assembly/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input.reads),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, key="out"),
        pairs=lambda wc, input, output: " ".join(",".join(x) for x in list(zip(input.reads, output.reads))),
        run_step=lambda wc: "t" if check_bool(wc, "Error_correction") else "f",
        prefilter=2,  # Ignore kmers with less than 2 occurance
        minprob=config["error_correction_minprob"],
        tossdepth=config["error_correction_minimum_kmer_depth"],
        tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
        lowdepthfraction=config["error_correction_lowdepth_fraction"],
        aggressive=config["error_correction_aggressive"],
        shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
    log:
        "{sample}/logs/assembly/pre_process/error_correction_{previous_steps}.log",
    benchmark:
        "logs/benchmarks/assembly/pre_process/{sample}_error_correction_{previous_steps}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ( 
        if [ '{params.run_step}' == 't' ]; then
            "tadpole.sh \
                prefilter={params.prefilter} \
                prealloc=1 \
                {params.inputs} \
                {params.outputs} \
                mode=correct \
                aggressive={params.aggressive} \
                tossjunk={params.tossjunk} \
                lowdepthfraction={params.lowdepthfraction} \
                tossdepth={params.tossdepth} \
                merge=t \
                shave={params.shave} \
                rinse={params.shave} \
                threads={threads} \
                pigz=t \
                unpigz=t \
                ecc=t \
                ecco=t \
                -Xmx{resources.java_mem}G \
        else
            echo 'NOTE: Skipping step and hard linking files instead.'
            mkdir -p "{output.reads}"
            touch "{output.stats}"
            for group in {params.pairs};
            do
                IFS=',' read -ra items <<< "$group"
                ln ${{items[0]}} ${{items[1]}}
            done
        fi
        ) 1>{log} 2>&1
        """


def get_assembly_command(sample_id, assembly_input):
    """
    Builds the assembly command using the file paths provided by the
    'run_assembly' rule's input directive.
    """
    assembler = sampleTable.loc[sample_id, 'Assembler']
    output_dir = f"results/assemblies/{sample_id}_{assembler}"

    # MEGAHIT / SPADES (Short Reads)
    if assembler == 'spades':
        cmd_part = "--meta"
        # Check for named inputs 'r1'/'r2' or 'se'
        if 'r1' in assembly_input and 'r2' in assembly_input:
            cmd_part += f" -1 {assembly_input.r1} -2 {assembly_input.r2}"
        elif 'r1' in assembly_input and not 'r2' in assembly_input:
            cmd_part += f" -s {assembly_input.r1}"
        else:
            raise ValueError(f"No trimmed short reads found for assembler '{assembler}' and sample '{sample_id}'.")
        
        # Add long reads for SPADES hybrid assembly
        if assembler == 'spades' and 'lr' in assembly_input:
            long_read_type = sampleTable.loc[sample_id].get('LongReadType', 'pacbio')
            cmd_part += f" --{long_read_type} {assembly_input.long_reads}"
            
        return f"spades.py {cmd_part} -o {output_dir}" if assembler == 'spades' else f"{assembler} {cmd_part} -o {output_dir}"

    elif assembler == 'megahit':
        pass
    # FLYE / METAMDBG (Long Reads)
    elif assembler in ['flye', 'metaMDBG']:
        if 'long_reads' not in assembly_input:
            raise ValueError(f"{assembler.capitalize()} requires long reads for sample '{sample_id}'.")
        
        long_reads = assembly_input.long_reads
        if assembler == 'flye':
            long_read_type = sampleTable.loc[sample_id].get('LongReadType', 'pacbio').replace('pacbio', 'pacbio-raw').replace('nanopore', 'nano-raw')
            return f"flye --{long_read_type} {long_reads} --out-dir {output_dir}"
        else: # metaMDBG
            return f"metaMDBG -i {long_reads} -o {output_dir}"
    else:
        raise ValueError(f"Unknown assembler '{assembler}'.")


rule run_assembly:
    input:
        # Use the helper function to get the correct mix of trimmed short reads and raw long reads
        unpack(get_assembly_inputs)
    output:
        contigs = "results/assemblies/{sample}_{assembler}/scaffolds.fasta"
    params:
        # The assembly command function now references the named inputs from this rule
        command = lambda wildcards, input: get_assembly_command(wildcards.sample, input)
    threads: 16
    log:
        "logs/assembly/{sample}_{assembler}.log"
    shell:
        """
        ({params.command} --threads {threads}) > {log} 2>&1
        """

localrules:
    rename_assembler_output,

rule rename_assembler_output:
    input:
        "{{sample}}/assembly/{sequences}.fasta".format(
        sequences="scaffolds" if config["spades_use_scaffolds"] else "contigs"
        ),
    output:
        temp("{sample}/assembly/{sample}_raw_contigs.fasta"),
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit sort -l -r -w 0 {input} > {output}"


rule rename_contigs:
    input:
        "{sample}/assembly/{sample}_raw_contigs.fasta",
    output:
        fasta="{sample}/assembly/{sample}_prefilter_contigs.fasta",
        mapping_table="{sample}/assembly/old2new_contig_names.tsv",
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    log:
        "{sample}/logs/assembly/post_process/rename_and_filter_size.log",
    params:
        minlength=config["minimum_contig_length"],
    conda:
        "../envs/fasta.yaml"
    script:
        "../scripts/rename_assembly.py"


if config["filter_contigs"]:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            query=get_quality_controlled_reads,
            target=rules.rename_contigs.output,
        output:
            bam=temp("{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam"),
        params:
            extra="-x sr",
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log",
        threads: config["simplejob_threads"]
        resources:
            mem_mb=config["simplejob_memory"] * 1000,
            time=config["simplejob_runtime"],
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

    rule pileup_prefilter:
        input:
            fasta="{sample}/assembly/{sample}_prefilter_contigs.fasta",
            bam="{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam",
        output:
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        params:
            pileup_secondary="t",
            minmapq=config["minimum_map_quality"],
        log:
            "{sample}/logs/assembly/post_process/pilup_prefilter_contigs.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["simplejob_threads"]
        resources:
            mem_mb=config["simplejob_memory"] * 1000,
            java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
            time=config["simplejob_runtime"],
        shell:
            "pileup.sh ref={input.fasta} in={input.bam} "
            " threads={threads} "
            " -Xmx{resources.java_mem}G "
            " covstats={output.covstats} "
            " concise=t "
            " minmapq={params.minmapq} "
            " secondary={params.pileup_secondary} "
            " 2> {log}"

    rule filter_by_coverage:
        input:
            fasta="{sample}/assembly/{sample}_prefilter_contigs.fasta",
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        output:
            fasta="{sample}/assembly/{sample}_final_contigs.fasta",
            removed_names="{sample}/assembly/{sample}_discarded_contigs.fasta",
        params:
            minc=config["minimum_average_coverage"],
            minp=config["minimum_percent_covered_bases"],
            minr=config.get("minimum_mapped_reads", MINIMUM_MAPPED_READS),
            minl=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
            trim=config.get("contig_trim_bp", CONTIG_TRIM_BP),
        log:
            "{sample}/logs/assembly/post_process/filter_by_coverage.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["simplejob_threads"]
        resources:
            mem=config["simplejob_memory"],
            java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
            time=config["simplejob_runtime"],
        shell:
            """filterbycoverage.sh in={input.fasta} \
            cov={input.covstats} \
            out={output.fasta} \
            outd={output.removed_names} \
            minc={params.minc} \
            minp={params.minp} \
            minr={params.minr} \
            minl={params.minl} \
            trim={params.trim} \
            -Xmx{resources.java_mem}G 2> {log}"""


# HACK: this makes two copies of the same file


else:  # no filter

    localrules:
        do_not_filter_contigs,

    rule do_not_filter_contigs:
        input:
            "{sample}/assembly/{sample}_prefilter_contigs.fasta",
        output:
            "{sample}/assembly/{sample}_final_contigs.fasta",
        shell:
            "cp {input} {output}"


localrules:
    finalize_contigs,


rule finalize_contigs:
    input:
        "{sample}/assembly/{sample}_final_contigs.fasta",
    output:
        "Assembly/fasta/{sample}.fasta",
    shell:
        "cp {input} {output}"


rule calculate_contigs_stats:
    input:
        get_assembly,
    output:
        "{sample}/assembly/contig_stats/final_contig_stats.txt",
    conda:
        "../envs/required_packages.yaml"
    log:
        "{sample}/logs/assembly/post_process/contig_stats_final.log",
    shell:
        "stats.sh in={input} format=3 out={output} &> {log}"


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_final_contigs:
    input:
        query=get_quality_controlled_reads,
        target="Assembly/fasta/{sample_contigs}.fasta",
    output:
        bam=temp("{sample_contigs}/sequence_alignment/{sample}.bam"),
    params:
        extra="-x sr",
        sorting="coordinate",
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/align_reads_to_filtered_contigs/{sample}_to_{sample_contigs}.txt"
    log:
        "{sample_contigs}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log",
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    wrapper:
        "v1.19.0/bio/minimap2/aligner"


rule pileup_contigs_sample:
    input:
        fasta=get_assembly,
        bam="{sample}/sequence_alignment/{sample}.bam",
    output:
        covhist="{sample}/assembly/contig_stats/postfilter_coverage_histogram.txt",
        covstats="{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
        bincov="{sample}/assembly/contig_stats/postfilter_coverage_binned.txt",
    params:
        pileup_secondary=(
            "t"
            if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS)
            else "f"
        ),
        minmapq=config["minimum_map_quality"],
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sample}.txt"
    log:
        "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",  # This log file is uesd for report
    conda:
        "../envs/required_packages.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
        time=config["simplejob_runtime"],
    shell:
        "pileup.sh "
        " ref={input.fasta} "
        " in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " hist={output.covhist} "
        " concise=t "
        " minmapq={params.minmapq} "
        " secondary={params.pileup_secondary} "
        " bincov={output.bincov} "
        " 2> {log} "


rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    conda:
        "../envs/required_packages.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        "samtools index {input}"


rule predict_genes:
    input:
        get_assembly,
    output:
        fna="{sample}/annotation/predicted_genes/{sample}.fna",
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
        gff="{sample}/annotation/predicted_genes/{sample}.gff",
    conda:
        "../envs/prodigal.yaml"
    log:
        "{sample}/logs/gene_annotation/prodigal.txt",
    benchmark:
        "logs/benchmarks/prodigal/{sample}.txt"
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """


localrules:
    get_contigs_from_gene_names,


rule get_contigs_from_gene_names:
    input:
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        tsv="{sample}/annotation/predicted_genes/{sample}.tsv",
    run:
        header = [
            "gene_id",
            "Contig",
            "Gene_nr",
            "Start",
            "Stop",
            "Strand",
            "Annotation",
        ]
        with open(output.tsv, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            with open(input.faa) as fin:
                gene_idx = 0
                for line in fin:
                    if line[0] == ">":
                        text = line[1:].strip().split(" # ")
                        old_gene_name = text[0]
                        text.remove(old_gene_name)
                        old_gene_name_split = old_gene_name.split("_")
                        gene_nr = old_gene_name_split[-1]
                        contig_nr = old_gene_name_split[-2]
                        sample = "_".join(
                            old_gene_name_split[: len(old_gene_name_split) - 2]
                        )
                        tsv.write(
                            "{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                                text="\t".join(text),
                                gene_id=old_gene_name,
                                i=gene_idx,
                                sample=sample,
                                gene_nr=gene_nr,
                                contig_nr=contig_nr,
                            )
                        )
                        gene_idx += 1



localrules:
    build_assembly_report,
    combine_contig_stats,


rule combine_contig_stats:
    input:
        contig_stats=expand(
            "{sample}/assembly/contig_stats/final_contig_stats.txt", sample=SAMPLES
        ),
        gene_tables=expand(
            "{sample}/annotation/predicted_genes/{sample}.tsv", sample=SAMPLES
        ),
        mapping_logs=expand(
            "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",
            sample=SAMPLES,
        ),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams=expand("{sample}/sequence_alignment/{sample}.bam", sample=SAMPLES),
    output:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/assembly/combine_contig_stats.log",
    script:
        "../scripts/combine_contig_stats.py"


rule build_assembly_report:
    input:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    output:
        report="reports/assembly_report.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/assembly/report.log",
    script:
        "../report/assembly_report.py"
