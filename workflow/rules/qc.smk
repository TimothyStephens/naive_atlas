import os
import re
import sys
from glob import glob
import warnings


def get_input_fastq(wildcards):
    """
    Get reads for QC by checking which files were provided.
    
    if sample has:
	R1       -> {reads:[R1]}
        R1,R2    -> {reads:[R1, R2]}
	R1,R2,LR -> {reads:[R1, R2], lr_reads:[LR]}
        LR       -> {lr_reads:[LR]}
    """
    sampleTable_info = sampleTable.loc[wildcards.sample, ].dropna()
    
    reads = {}
    headers = []
    if 'Reads_raw_R1' in sampleTable_info:
        headers.append('Reads_raw_R1')
    if 'Reads_raw_R2' in sampleTable_info:
        headers.append('Reads_raw_R2')
    if headers:
        reads["reads"] = get_files_from_sampleTable(wildcards.sample, headers)
    
    if 'Reads_raw_Long' in sampleTable_info: # If only LR provided
        headers.append('Reads_raw_Long')
        reads["lr_reads"] = get_files_from_sampleTable(wildcards.sample, ['Reads_raw_Long'])
    
    if not headers:
        ValueError(f"No reads found for sample '{wildcards.sample}'.")
    
    return(reads)

def get_output_fastq(wildcards, step):
    #print(f"wildcards: {wildcards}; step: {step}")
    files = expand(
        "samples/{{sample}}/sequence_quality_control/{step}/{{sample}}_{fraction}.fastq.gz",
        step=step,
        fraction=get_fractions(wildcards.sample),
    )
    return(files)

def check_interleaved(sample):
    return(sampleTable.loc[sample, "Interleaved"])

def check_paired(sample):
    return(
        len(sampleTable.loc[sample, ["Reads_raw_R1", "Reads_raw_R2"]].dropna()) == 2 or 
        (len(sampleTable.loc[sample, ["Reads_raw_R1"]].dropna()) == 1 and sampleTable.loc[sample, "Interleaved"])
    )

def check_bool(sample, column):
    try:
        b = sampleTable.loc[sample, column].item()
    except KeyError as e:
        raise KeyError(f"Sample '{sample}' or column '{column}' are missing from sampleTable. {e}")
    
    if not isinstance(b, bool):
        raise KeyError(f"Column '{column}' did not return a bool value, it returned '{b}'.")
    
    return(b)

def get_fractions(sample):
    sampleTable_info = sampleTable.loc[sample, ].dropna()
    
    fractions = []
    if 'Reads_raw_R1' in sampleTable_info and 'Reads_raw_R2' in sampleTable_info: 
        fractions.extend(["R1", "R2"])
    elif 'Reads_raw_R1' in sampleTable_info and not 'Reads_raw_R2' in sampleTable_info:
        if sampleTable.loc[sample, "Interleaved"]:
            fractions.extend(["R1", "R2"])
        else:
            fractions.extend(["SE"])
    if 'Reads_raw_Long' in sampleTable_info:
        fractions.extend(["LR"])
    
    return(fractions)


# List QC steps performed. Used to easily run stats once we are done.
PROCESSED_STEPS = []



####
#### Raw
####
PROCESSED_STEPS.append("1_raw")

rule initialize_qc_PE:
    input:
        unpack(get_input_fastq),
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/1_raw_reads_R1.fastq.gz", 
            "samples/{sample}/sequence_quality_control/cleaning/1_raw_reads_R2.fastq.gz"
        ]),
    priority: 80
    params:
        inputs =lambda wc, input:  io_params_for_tadpole(input.reads,  "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        interleaved=lambda wc: "t" if check_interleaved(wc.sample) else "f",
        verifypaired="t",
        extra=config["importqc_params"],
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/1_raw_PE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/1_raw_PE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "account"),
    shell:
        """
        (
        reformat.sh \\
            {params.inputs} \\
            interleaved={params.interleaved} \\
            {params.outputs} \\
            {params.extra} \\
            overwrite=true \\
            verifypaired={params.verifypaired} \\
            threads={threads} \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """


rule initialize_qc_SE:
    input:
        unpack(get_input_fastq),
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/1_raw_reads_SE.fastq.gz",
        ]),
    priority: 80
    params:
        inputs =lambda wc, input:  io_params_for_tadpole(input.reads,  "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        interleaved="f",
        verifypaired="f",
        extra=config["importqc_params"],
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/1_raw_SE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/1_raw_SE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "account"),
    shell:
        """
        (
        reformat.sh \\
            {params.inputs} \\
            interleaved={params.interleaved} \\
            {params.outputs} \\
            {params.extra} \\
            overwrite=true \\
            verifypaired={params.verifypaired} \\
            threads={threads} \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """


rule initialize_qc_LR:
    input:
        unpack(get_input_fastq),
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/1_raw_reads_LR.fastq.gz",
        ]),
    priority: 80
    params:
        inputs =lambda wc, input:  io_params_for_tadpole(input.lr_reads,  "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        interleaved="f",
        verifypaired="f",
        extra=config["importqc_params"],
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/1_raw_LR.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/1_raw_LR.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "initialize_qc", "account"),
    shell:
        """
        (
        reformat.sh \\
            {params.inputs} \\
            interleaved={params.interleaved} \\
            {params.outputs} \\
            {params.extra} \\
            overwrite=true \\
            verifypaired={params.verifypaired} \\
            threads={threads} \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """



####
#### DeDuplicated
####
PROCESSED_STEPS.append("2_deduplicated")

def deduplicate_reads_command(inputs, outputs, outdir, pairs, run_step, dupesubs, only_optical, threads, resources):
    if run_step == 't':
        cmd = f"""
        clumpify.sh \\
            {inputs} \\
            {outputs} \\
            overwrite=true \\
            dedupe=t \\
            dupesubs={dupesubs} \\
            optical={only_optical} \\
            threads={threads} \\
            pigz=t unpigz=t \\
            -Xmx{resources.java_mem}M
        """
    else:
        cmd = f"""
        echo 'Skipping step, hard linking files instead.'
        mkdir -p "{outdir}"
        IFS=';' read -ra groups <<< '{pairs}'
        for group in "${{groups[@]}}";
        do
            IFS=',' read -ra items <<< "$group"
            cp "${{items[0]}}" "${{items[1]}}"
        done
        """
    return(cmd)


rule deduplicate_reads_PE:
    input:
        reads=rules.initialize_qc_PE.output.reads,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_reads_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_reads_R2.fastq.gz"
        ]),
    params:
        command = lambda wc, input, output, threads, resources: deduplicate_reads_command(
            inputs =io_params_for_tadpole(input.reads,   "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            outdir=f"samples/{wc.sample}/sequence_quality_control/2_deduplicated",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "DeDuplicate") else "f",
            dupesubs=config["duplicates_allow_substitutions"],
            only_optical=("t" if config["duplicates_only_optical"] else "f"),
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_PE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_PE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "deduplicate_reads", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule deduplicate_reads_SE:
    input:
        reads=rules.initialize_qc_SE.output.reads,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_reads_SE.fastq.gz",
        ]),
    params:
        command = lambda wc, input, output, threads, resources: deduplicate_reads_command(
            inputs =io_params_for_tadpole(input.reads,   "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            outdir=f"samples/{wc.sample}/sequence_quality_control/2_deduplicated",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "DeDuplicate") else "f",
            dupesubs=config["duplicates_allow_substitutions"],
            only_optical=("t" if config["duplicates_only_optical"] else "f"),
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_SE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_SE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "deduplicate_reads", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule deduplicate_reads_LR:
    input:
        reads=rules.initialize_qc_LR.output.reads,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_reads_LR.fastq.gz",
        ]),
    params:
        command = lambda wc, input, output, threads, resources: deduplicate_reads_command(
            inputs =io_params_for_tadpole(input.reads,   "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            outdir=f"samples/{wc.sample}/sequence_quality_control/2_deduplicated",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "DeDuplicate") else "f",
            dupesubs=config["duplicates_allow_substitutions"],
            only_optical=("t" if config["duplicates_only_optical"] else "f"),
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_LR.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/2_deduplicated_LR.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "deduplicate_reads", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "deduplicate_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """



####
#### Filtered
####
PROCESSED_STEPS.append("3_quality_filtered")


def apply_quality_filter_command(inputs, outputs, stats, outdir, pairs, run_step,
                                ref, mink, ktrim, trimq,hdist, k, qtrim, 
                                error_correction_pe, minlength, minbasefrequency, 
                                interleaved, maxns, prealloc, threads, resources):
    if run_step == 't':
        cmd = f"""
        bbduk.sh \\
            {inputs} \\
            {ref} \\
            interleaved={interleaved} \\
            {outputs} \\
            stats={stats} \\
            overwrite=true \\
            qout=33 \\
            trd=t \\
            {hdist} \\
            {k} \\
            {ktrim} \\
            {mink} \\
            trimq={trimq} \\
            qtrim={qtrim} \\
            threads={threads} \\
            minlength={minlength} \\
            maxns={maxns} \\
            minbasefrequency={minbasefrequency} \\
            ecco={error_correction_pe} \\
            prealloc={prealloc} \\
            pigz=t unpigz=t \\
            -Xmx{resources.java_mem}M
        """
    else:
        cmd = f"""
        echo 'Skipping step, hard linking files instead.'
        mkdir -p "{outdir}"
        touch "{stats}"
        IFS=';' read -ra groups <<< '{pairs}'
        for group in "${{groups[@]}}";
        do
            IFS=',' read -ra items <<< "$group"
            cp "${{items[0]}}" "${{items[1]}}"
        done
        """
    return(cmd)


rule apply_quality_filter_PE:
    input:
        reads=rules.deduplicate_reads_PE.output.reads,
        adapters=ADAPTERS,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_R2.fastq.gz"
        ]),
        stats="samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_PE_stats.txt",
    params:
        command = lambda wc, input, output, threads, resources: apply_quality_filter_command(
            inputs =io_params_for_tadpole(input.reads,  "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            stats=output.stats,
            outdir=f"samples/{wc.sample}/sequence_quality_control/3_quality_filtered",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Quality_filter") else "f",
            ref=(
                "ref=%s" % ADAPTERS
                if (ADAPTERS is not None)
                else ""
            ),
            mink="mink=%d" % config["preprocess_adapter_min_k"],
            ktrim="ktrim=%s" % config["preprocess_kmer_trim"],
            trimq=config["preprocess_minimum_base_quality"],
            hdist="hdist=%d" % config["preprocess_allowable_kmer_mismatches"],
            k="k=%d" % config["preprocess_reference_kmer_match_length"],
            qtrim=config["preprocess_qtrim"],
            error_correction_pe=(
                "t"
                if config["error_correction_overlapping_pairs"]
                else "f"
            ),
            minlength=config["preprocess_minimum_passing_read_length"],
            minbasefrequency=config["preprocess_minimum_base_frequency"],
            # we require the user to reformat to R1 and R2, non-interleaved files
            interleaved="f",
            maxns=config["preprocess_max_ns"],
            prealloc=config["preallocate_ram"],
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_PE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_PE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "apply_quality_filter", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule apply_quality_filter_SE:
    input:
        reads=rules.deduplicate_reads_SE.output.reads,
        adapters=ADAPTERS,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_SE.fastq.gz",
        ]),
        stats="samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_SE_stats.txt",
    params:
        command = lambda wc, input, output, threads, resources: apply_quality_filter_command(
            inputs =io_params_for_tadpole(input.reads,  "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            stats=output.stats,
            outdir=f"samples/{wc.sample}/sequence_quality_control/3_quality_filtered",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Quality_filter") else "f",
            ref=(
                "ref=%s" % ADAPTERS
                if (ADAPTERS is not None)
                else ""
            ),
            mink="mink=%d" % config["preprocess_adapter_min_k"],
            ktrim="ktrim=%s" % config["preprocess_kmer_trim"],
            trimq=config["preprocess_minimum_base_quality"],
            hdist="hdist=%d" % config["preprocess_allowable_kmer_mismatches"],
            k="k=%d" % config["preprocess_reference_kmer_match_length"],
            qtrim=config["preprocess_qtrim"],
            error_correction_pe="f",
            minlength=config["preprocess_minimum_passing_read_length"],
            minbasefrequency=config["preprocess_minimum_base_frequency"],
            # we require the user to reformat to R1 and R2, non-interleaved files
            interleaved="f",
            maxns=config["preprocess_max_ns"],
            prealloc=config["preallocate_ram"],
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_SE.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_SE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "apply_quality_filter", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule apply_quality_filter_LR:
    input:
        reads=rules.deduplicate_reads_LR.output.reads,
        adapters=ADAPTERS,
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_LR.fastq.gz",
        ]),
        stats="samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_reads_LR_stats.txt",
    params:
        command = lambda wc, input, output, threads, resources: apply_quality_filter_command(
            inputs =io_params_for_tadpole(input.reads,  "in"),
            outputs=io_params_for_tadpole(output.reads, "out"),
            stats=output.stats,
            outdir=f"samples/{wc.sample}/sequence_quality_control/3_quality_filtered",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Quality_filter") else "f",
            ref=(
                "ref=%s" % ADAPTERS
                if (ADAPTERS is not None)
                else ""
            ),
            mink="mink=%d" % config["preprocess_adapter_min_k"],
            ktrim="ktrim=%s" % config["preprocess_kmer_trim"],
            trimq=config["preprocess_minimum_base_quality"],  
            hdist="hdist=%d" % config["preprocess_allowable_kmer_mismatches"],
            k="k=%d" % config["preprocess_reference_kmer_match_length"],
            qtrim=config["preprocess_qtrim"],
            error_correction_pe="f",
            minlength=config["preprocess_minimum_passing_read_length"],
            minbasefrequency=config["preprocess_minimum_base_frequency"],
            # we require the user to reformat to R1 and R2, non-interleaved files
            interleaved="f",
            maxns=config["preprocess_max_ns"],
            prealloc=config["preallocate_ram"],
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_LR.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/cleaning/3_quality_filtered_LR.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "apply_quality_filter", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "apply_quality_filter", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """



####
#### Contaminant References
####
# if there are no references, decontamination will be skipped
if len(CONTAMINANT_REFERENCES.keys()) > 0:
    PROCESSED_STEPS.append("4_decontaminated")

    rule build_decontamination_db:
        input:
            ancient(CONTAMINANT_REFERENCES.values()),
        output:
            "ref/genome/1/summary.txt",
        params:
            k=config["contaminant_kmer_length"],
            refs_in=" ".join(
                [
                    "ref_%s=%s" % (n, fa)
                    for n, fa in CONTAMINANT_REFERENCES.items()
                ]
            ),
        log:
            "logs/qc/sequence_quality_control/cleaning/4_decontamination_build_db.log",
        benchmark:
            "benchmarks/qc/sequence_quality_control/cleaning/4_decontamination_build_db.tsv",
        conda:
            "../envs/required_packages.yaml"
        threads: lambda wc: get_resource(wc, None, 1, "build_decontamination_db", "threads")
        resources:
            mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "build_decontamination_db", "mem_mb"),
            java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "build_decontamination_db", "java_mem"),
            runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "build_decontamination_db", "time_min"),
            slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "build_decontamination_db", "partition"),
            slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "build_decontamination_db", "account"),
        shell:
            """
            (
            bbsplit.sh \\
                -Xmx{resources.java_mem}M \\
                {params.refs_in} \\
                threads={threads} \\
                k={params.k} \\
                local=t
            ) 1>{log} 2>&1
            """
    
    
    def run_decontamination_command(inputs, outputs, stats, contaminant_folder, outdir, pairs, run_step, pacbio,
                                    paired, maxindel, minratio, minhits, ambiguous, k, threads, resources):
        if run_step == 't' and not pacbio:
            cmd = f"""
            bbsplit.sh \\
                {inputs} \\
                {outputs} \\
                basename={contaminant_folder}/%_R#.fastq.gz \\
                maxindel={maxindel} \\
                minratio={minratio} \\
                minhits={minhits} \\
                ambiguous={ambiguous} \\
                refstats={stats} \\
                threads={threads} \\
                k={k} \\
                local=t \\
                machineout=t \\
                pigz=t unpigz=t ziplevel=9 \\
                -Xmx{resources.java_mem}M
            """
        elif run_step == 't' and pacbio:
            cmd = f"""
            reformat.sh \\
                {inputs} \\
                out=stdout.fa \\
            | bbsplit.sh \\
                in=stdin.fa \\
                basename={contaminant_folder}/%_LR.fastq.gz \\
                maxindel={maxindel} \\
                minratio={minratio} \\
                minhits={minhits} \\
                ambiguous={ambiguous} \\
                refstats={stats} \\
                threads={threads} \\
                k={k} \\
                local=t \\
                machineout=t \\
                pigz=t unpigz=t ziplevel=9 \\
                -Xmx{resources.java_mem}M \\
            && zcat {contaminant_folder}/*_LR.fastq.gz \\
              | reformat.sh in=stdin.fq out=stdout.fa int=f minlength=50 \\
              | grep '>' \\
              | sed -e 's/>//' -e 's@_part_[0-9+]\+$@@' \\
              | sort | uniq \\
              > {contaminant_folder}/contaminant_read_names.txt \\
            && filterbyname.sh \\
                names={contaminant_folder}/contaminant_read_names.txt \\
                {inputs} \\
                {outputs}
            """
        else:
            cmd = f"""
            echo 'Skipping step, hard linking files instead.'
            mkdir -p "{outdir}"
            touch "{stats}"
            IFS=';' read -ra groups <<< '{pairs}'
            for group in "${{groups[@]}}";
            do
                IFS=',' read -ra items <<< "$group"
                cp "${{items[0]}}" "${{items[1]}}"
            done
            """
        return(cmd)
    
    
    rule run_decontamination_PE:
        input:
            reads=rules.apply_quality_filter_PE.output.reads,
            db="ref/genome/1/summary.txt",
        output:
            reads=temp([
                "samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_R1.fastq.gz",
                "samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_R2.fastq.gz"
            ]),
            stats="samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_PE_stats.txt",
            contaminant_folder=directory("samples/{sample}/sequence_quality_control/cleaning/4_contaminated_reads_PE")
        params:
            command = lambda wc, input, output, threads, resources: run_decontamination_command(
                inputs =io_params_for_tadpole(input.reads,  "in"),
                outputs=io_params_for_tadpole(output.reads, "outu"),
                stats=output.stats,
                contaminant_folder=output.contaminant_folder,
                outdir=output.contaminant_folder,
                pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
                run_step="t" if check_bool(wc, "Remove_contaminants") else "f",
                pacbio=False,
                paired="true",
                maxindel=config["contaminant_max_indel"],
                minratio=config["contaminant_min_ratio"],
                minhits=config["contaminant_minimum_hits"],
                ambiguous=config["contaminant_ambiguous"],
                k=config["contaminant_kmer_length"],
                threads=threads,
                resources=resources
            )
        log:
            "logs/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_PE.log",
        benchmark:
            "benchmarks/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_PE.tsv",
        conda:
            "../envs/required_packages.yaml"
        threads: lambda wc: get_resource(wc, None, 1, "run_decontamination", "threads")
        resources:
            mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "mem_mb"),
            java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "java_mem"),
            runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "time_min"),
            slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "partition"),
            slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "account"),
        shell:
            """
            ({params.command}) 1>{log} 2>&1
            """


    rule run_decontamination_SE:
        input:
            reads=rules.apply_quality_filter_SE.output.reads,
            db="ref/genome/1/summary.txt",
        output:
            reads=temp([
                "samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_SE.fastq.gz",
            ]),
            stats="samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reference_SE_stats.txt",
            contaminant_folder=directory("samples/{sample}/sequence_quality_control/cleaning/4_contaminated_reads_SE")
        params:
            command = lambda wc, input, output, threads, resources: run_decontamination_command(
                inputs =io_params_for_tadpole(input.reads,  "in"),
                outputs=io_params_for_tadpole(output.reads, "outu"),
                stats=output.stats,
                contaminant_folder=output.contaminant_folder,
                outdir=output.contaminant_folder,
                pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
                run_step="t" if check_bool(wc, "Remove_contaminants") else "f",
                pacbio=False,
                paired="false",
                maxindel=config["contaminant_max_indel"],
                minratio=config["contaminant_min_ratio"],
                minhits=config["contaminant_minimum_hits"],
                ambiguous=config["contaminant_ambiguous"],
                k=config["contaminant_kmer_length"],
                threads=threads,
                resources=resources
            )
        log:
            "logs/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_SE.log",
        benchmark:
            "benchmarks/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_SE.tsv",
        conda:
            "../envs/required_packages.yaml"
        threads: lambda wc: get_resource(wc, None, 1, "run_decontamination", "threads")
        resources:
            mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "mem_mb"),
            java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "java_mem"),
            runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "time_min"),
            slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "partition"),
            slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "account"),
        shell:
            """
            ({params.command}) 1>{log} 2>&1
            """


    rule run_decontamination_LR:
        input:
            reads=rules.apply_quality_filter_LR.output.reads,
            db="ref/genome/1/summary.txt",
        output:
            reads=temp([
                "samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_LR.fastq.gz",
            ]),
            stats="samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_reads_LR_stats.txt",
            contaminant_folder=directory("samples/{sample}/sequence_quality_control/cleaning/4_contaminated_reads_LR")
        params:
            command = lambda wc, input, output, threads, resources: run_decontamination_command(
                inputs =io_params_for_tadpole(input.reads,  "in"),
                outputs=io_params_for_tadpole(output.reads, "out"),
                stats=output.stats,
                contaminant_folder=output.contaminant_folder,
                outdir=output.contaminant_folder,
                pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
                run_step="t" if check_bool(wc, "Remove_contaminants") else "f",
                pacbio=True,
                #run_step="f",
                paired="false",
                maxindel=config["contaminant_max_indel"],
                minratio=config["contaminant_min_ratio"],
                minhits=config["contaminant_minimum_hits"],
                ambiguous=config["contaminant_ambiguous"],
                k=config["contaminant_kmer_length"],
                threads=threads,
                resources=resources
            )
        log:
            "logs/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_LR.log",
        benchmark:
            "benchmarks/samples/{sample}/sequence_quality_control/cleaning/4_decontaminated_LR.tsv",
        conda:
            "../envs/required_packages.yaml"
        threads: lambda wc: get_resource(wc, None, 1, "run_decontamination", "threads")
        resources:
            mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "mem_mb"),
            java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "java_mem"),
            runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "time_min"),
            slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "partition"),
            slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_decontamination", "account"),
        shell:
            """
            ({params.command}) 1>{log} 2>&1
            """



####
#### QC
####
PROCESSED_STEPS.append("5_final")

localrules:
    qcreads_PE,
    qcreads_SE,
    qcreads_LR

rule qcreads_PE:
    input:
        reads=(
            rules.run_decontamination_PE.output.reads
            if len(CONTAMINANT_REFERENCES.keys()) > 0
            else rules.apply_quality_filter_PE.output.reads
        ),
    output:
        reads=temp([   
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_R2.fastq.gz"
        ]),
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/5_qcreads_PE.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil
        for i, f in enumerate(input.reads):
            shutil.copy(f, output.reads[i])


rule qcreads_SE:
    input:
        reads=(   
            rules.run_decontamination_SE.output.reads
            if len(CONTAMINANT_REFERENCES.keys()) > 0
            else rules.apply_quality_filter_SE.output.reads
        ),
    output:
        reads=temp([
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_SE.fastq.gz",
        ]),
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/5_qcreads_SE.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil
        for i, f in enumerate(input.reads):
            shutil.copy(f, output.reads[i])


rule qcreads_LR:
    input:
        reads=(
            rules.run_decontamination_LR.output.reads
            if len(CONTAMINANT_REFERENCES.keys()) > 0
            else rules.apply_quality_filter_LR.output.reads
        ),
    output:
        reads=temp([   
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_LR.fastq.gz",
        ]),
    log:
        "logs/samples/{sample}/sequence_quality_control/cleaning/5_qcreads_LR.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil
        for i, f in enumerate(input.reads):
            shutil.copy(f, output.reads[i])



####
#### Copy Cleaned Reads
####
localrules:
    copy_reads_PE,
    copy_reads_SE,
    copy_reads_LR


rule copy_reads_PE:
    input:
        reads=[
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_R2.fastq.gz"
        ],
    output:
        reads=[
            "samples/{sample}/sequence_quality_control/{sample}_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/{sample}_R2.fastq.gz"
        ],
    log:
        "logs/samples/{sample}/sequence_quality_control/5_copy_reads_PE.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil, os
        import pandas as pd
        
        for i in range(len(output.reads)):
            with open(output.reads[i], "wb") as outFile, open(input.reads[i], "rb") as infile:
                shutil.copyfileobj(infile, outFile)


rule copy_reads_SE:
    input:
        reads=[
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_SE.fastq.gz",
        ],
    output:
        reads=[
            "samples/{sample}/sequence_quality_control/{sample}_SE.fastq.gz",
        ],
    log:
        "logs/samples/{sample}/sequence_quality_control/5_copy_reads_SE.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil, os
        import pandas as pd
        
        for i in range(len(output.reads)):
            with open(output.reads[i], "wb") as outFile, open(input.reads[i], "rb") as infile:
                shutil.copyfileobj(infile, outFile)


rule copy_reads_LR:
    input:
        reads=[
            "samples/{sample}/sequence_quality_control/cleaning/5_final_reads_LR.fastq.gz",
        ],
    output:
        reads=[
            "samples/{sample}/sequence_quality_control/{sample}_LR.fastq.gz",
        ],
    log:
        "logs/samples/{sample}/sequence_quality_control/5_copy_reads_LR.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import shutil, os
        import pandas as pd
        
        for i in range(len(output.reads)):
            with open(output.reads[i], "wb") as outFile, open(input.reads[i], "rb") as infile:
                shutil.copyfileobj(infile, outFile)


def get_quality_controlled_reads(wildcards, as_dict=False, subset=None):
    if subset is None:
        prefix=f"samples/{wildcards.sample}/sequence_quality_control/{wildcards.sample}"
    else:
        prefix=f"samples/{wildcards.sample}/sequence_quality_control/cleaning/{subset}_reads"
    if as_dict:
        files = {}
        for fraction in get_fractions(wildcards.sample):
            files[fraction] = f"{prefix}_{fraction}.fastq.gz"
    else:
        files = []
        for fraction in get_fractions(wildcards.sample):
            files.append(f"{prefix}_{fraction}.fastq.gz")
    return(files)



####
#### STATS
####

#
# Read counts
#
rule get_read_counts:
    input:
        unpack(lambda wc: get_quality_controlled_reads(wc, as_dict=True, subset=f"{wc.step}")),
    output:
        zipped_dir="samples/{sample}/sequence_quality_control/read_stats/{step}.zip",
        read_counts=temp("samples/{sample}/sequence_quality_control/read_stats/{step}_read_counts.tsv"),
    params:
        folder=lambda wc, output: os.path.splitext(output['zipped_dir'])[0],
    log:
        "logs/samples/{sample}/sequence_quality_control/read_stats/{step}.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/read_stats/{step}.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "get_read_counts", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_counts", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_counts", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_counts", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_counts", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_counts", "account"),
    priority: 30
    script:
        "../scripts/get_read_counts.py"


localrules:
    write_read_counts,
    combine_read_counts,


rule write_read_counts:
    input:
        read_count_files=expand(
            "samples/{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv",
            step=PROCESSED_STEPS,
        ),
    output:
        read_stats="samples/{sample}/sequence_quality_control/read_stats/read_counts.tsv",
    log:
        "logs/samples/{sample}/sequence_quality_control/read_stats/write_read_counts.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        from utils.io import pandas_concat
        
        pandas_concat(
            list(input.read_count_files),
            output.read_stats,
            sep="\t",
            index_col=[0, 1],
            axis=0,
        )


rule combine_read_counts:
    input:
        expand(
            "samples/{sample}/sequence_quality_control/read_stats/read_counts.tsv",
            sample=SAMPLES,
        ),
    output:
        "stats/read_counts.tsv",
    log:
        "logs/qc/sequence_quality_control/read_stats/combine_read_counts.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        from utils.io import pandas_concat
        
        pandas_concat(list(input), output[0], sep="\t", index_col=[0, 1], axis=0)



#
# Read length and insert size
#
rule get_read_length_hist:
    input:
        unpack(lambda wc: get_quality_controlled_reads(wc, as_dict=True)),
    output:
        insertHist_pe="samples/{sample}/sequence_quality_control/read_stats/QC_insert_size_hist_PE.txt",
        insertHist_se="samples/{sample}/sequence_quality_control/read_stats/QC_insert_size_hist_SE.txt",
        insertHist_lr="samples/{sample}/sequence_quality_control/read_stats/QC_insert_size_hist_LR.txt",
        lenHist_pe="samples/{sample}/sequence_quality_control/read_stats/QC_read_length_hist_PE.txt",
        lenHist_se="samples/{sample}/sequence_quality_control/read_stats/QC_read_length_hist_SE.txt",
        lenHist_lr="samples/{sample}/sequence_quality_control/read_stats/QC_read_length_hist_LR.txt",
    params:
        kmer=config["merging_k"],
        extend2=config["merging_extend2"],
        flags="loose ecct",
        minprob=config["bbmerge_minprob"],
    log:
        "logs/samples/{sample}/sequence_quality_control/read_stats/calculate_read_length.log",
    benchmark:
        "benchmarks/samples/{sample}/sequence_quality_control/read_stats/calculate_read_length.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "get_read_length_hist", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_length_hist", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_length_hist", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_length_hist", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_length_hist", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "get_read_length_hist", "account"),
    script:
        "../scripts/get_read_length.py"



localrules:
    combine_read_length_hist,
    combine_insert_hist,


rule combine_read_length_hist:
    input:
        expand(
            "samples/{sample}/sequence_quality_control/read_stats/QC_read_length_hist_{fraction}.txt",
            sample=SAMPLES,
            fraction=["PE", "SE", "LR"]
        ),
    output:
        "stats/read_length_stats.tsv",
    log:
        "logs/qc/sequence_quality_control/read_stats/combine_read_length_hist.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()
        
        i = 0
        for length_file in input:
            i+=1
            sample = length_file.split(os.path.sep)[1]
            fraction = length_file.split(os.path.sep)[-1].replace("QC_read_length_hist_", "").replace(".txt", "")
            data = parse_comments(length_file)
            data = pd.Series(data)[
                ["Reads", "Bases", "Max", "Min", "Avg", "Median", "Mode", "Std_Dev"]
            ]
            data = pd.concat([pd.Series([sample, fraction], index=["sample", "fraction"]), data])
            stats[str(i)] = data

        stats.to_csv(output[0], sep="\t")


rule combine_insert_hist:
    input:
        expand(
            "samples/{sample}/sequence_quality_control/read_stats/QC_insert_size_hist_{fraction}.txt",
            sample=SAMPLES,
            fraction=["PE", "SE", "LR"]
        ),
    output:
        "stats/insert_stats.tsv",
    log:
        "logs/qc/sequence_quality_control/read_stats/combine_insert_hist.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()
        
        i = 0
        for insert_file in input:
            i+=1
            sample = insert_file.split(os.path.sep)[1]
            fraction = insert_file.split(os.path.sep)[-1].replace("QC_insert_size_hist_", "").replace(".txt", "")
            data = parse_comments(insert_file)
            data = pd.Series(data)[
                ["Mean", "Median", "Mode", "STDev", "PercentOfPairs"]
            ]
            data = pd.concat([pd.Series([sample, fraction], index=["sample", "fraction"]), data])
            stats[str(i)] = data

        stats.T.to_csv(output[0], sep="\t")



####
#### Build QC Report
####

localrules:
    build_qc_report

rule build_qc_report:
    input:
        zipfiles_QC=expand(
            "samples/{sample}/sequence_quality_control/read_stats/5_final.zip",
            sample=SAMPLES
        ),
        read_counts="stats/read_counts.tsv",
        read_length_stats="stats/read_length_stats.tsv",
        read_insert_size_stats="stats/insert_stats.tsv",
    output:
        report="reports/qc_report.html",
    log:
        "logs/qc/sequence_quality_control/read_stats/report.log",
    params:
        min_quality=config["preprocess_minimum_base_quality"],
        samples=SAMPLES,
    conda:
        "../envs/report.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    script:
        "../report/qc_report.py"


####
#### Done
####

localrules:
    finalize_sample_qc

rule finalize_sample_qc:
    input:
        report="reports/qc_report.html",
    output:
        flag=touch("samples/{sample}/sequence_quality_control/finished_QC"),
    log:
        "logs/samples/{sample}/sequence_quality_control/finalize_sample_qc.log",
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),

