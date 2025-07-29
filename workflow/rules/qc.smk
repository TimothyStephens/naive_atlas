import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules:
    build_qc_report,
    combine_read_length_stats,
    combine_read_counts,


def get_ribosomal_rna_input(wildcards):
    data_type = config["data_type"]

    clean_reads = get_output_fastq(wildcards, "4_decontamination")
    rrna_reads = expand(
        "{sample}/sequence_quality_control/4_contaminants/rRNA_{fraction}.fastq.gz",
        fraction=get_fractions(wildcards.sample),
        sample=wildcards.sample,
    )

    if data_type == "metagenome" and "rrna" in [
        c.lower() for c in config["contaminant_references"].keys()
    ]:
        return {"clean_reads": clean_reads, "rrna_reads": rrna_reads}
    else:
        return {"clean_reads": clean_reads}


def get_input_fastq(wildcards, return_lr=False):
    """
    Get reads for QC by checking which files were provided.
    
    if sample has:
        R1       -> R1
        R1,R2    -> R1,R2
        R1,R2,LR -> R1,R2 (assume LR are for scaffolding only)
        LR       -> LR (assume LR is high quality or coverage for LR-only assembly)
    """
    sampleTable_info = sampleTable.loc[wildcards.sample, ].dropna()
    
    headers = []
    if not return_lr:
        if 'Reads_raw_R1' in sampleTable_info:
            headers.append('Reads_raw_R1')
        if 'Reads_raw_R2' in sampleTable_info:
            headers.append('Reads_raw_R2')
    else:
        if 'Reads_raw_Long' in sampleTable_info: # If only LR provided
            headers.append('Reads_raw_Long')
    
    if not headers:
        ValueError(f"No reads found for sample '{wildcards.sample}'.")
    
    return get_files_from_sampleTable(wildcards.sample, headers)

def get_output_fastq(wildcards, step):
    #print(f"wildcards: {wildcards}; step: {step}")
    files = expand(
        "{sample}/sequence_quality_control/{step}/{sample}_{fraction}.fastq.gz",
        sample=wildcards.sample,
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
    
    if 'Reads_raw_R1' in sampleTable_info and 'Reads_raw_R2' in sampleTable_info: 
        fractions = ["R1", "R2"]
    elif 'Reads_raw_R1' in sampleTable_info and not 'Reads_raw_R2' in sampleTable_info:
        if sampleTable.loc[sample, "Interleaved"]:
            fractions = ["R1", "R2"]
        else:
            fractions = ["SE"]
    elif not 'Reads_raw_R1' in sampleTable_info and not 'Reads_raw_R2' in sampleTable_info and 'Reads_raw_Long' in sampleTable_info:
        fractions = ["LR"]
    
    return(fractions)


# List QC steps performed. Used to easily run stats once we are done.
PROCESSED_STEPS = []


####
#### Raw
####
PROCESSED_STEPS.append("1_raw")

rule initialize_qc:
    input:
        reads=get_input_fastq,
    output:
        reads=temp(lambda wc: get_output_fastq(wc, "1_raw")),
    priority: 80
    params:
        inputs =lambda wc, input:  io_params_for_tadpole(input.reads,  "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        interleaved=lambda wc: "t" if check_interleaved(wc.sample) else "f",
        verifypaired=lambda wc: "t" if check_paired(wc.sample) or check_interleaved(wc.sample) else "f",
        extra=config["importqc_params"],
    log:
        "{sample}/logs/QC/1_raw.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["simplejob_threads"]
    benchmark:
        "logs/benchmarks/QC/1_raw/{sample}.txt"
    resources:
        mem=config["simplejob_memory"],
        java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
    shell:
        "reformat.sh "
        " {params.inputs} "
        " interleaved={params.interleaved} "
        " {params.outputs} "
        " {params.extra} "
        " overwrite=true "
        " verifypaired={params.verifypaired} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " 1>{log} 2>&1 "


####
#### DeDuplicated
####
PROCESSED_STEPS.append("2_deduplicated")

rule deduplicate_reads:
    input:
        reads=lambda wc: get_output_fastq(wc, "1_raw"),
    output:
        reads=temp(lambda wc: get_output_fastq(wc, "2_deduplicated")),
    params:
        inputs =lambda wc, input: io_params_for_tadpole(input.reads,   "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        pairs=lambda wc, input, output: " ".join(",".join(x) for x in list(zip(input.reads, output.reads))),
        run_step=lambda wc: "t" if check_bool(wc, "DeDuplicate") else "f",
        dupesubs=config["duplicates_allow_substitutions"],
        only_optical=("t" if config.get("duplicates_only_optical") else "f"),
    log:
        "{sample}/logs/QC/2_deduplicated.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["simplejob_threads"]
    benchmark:
        "logs/benchmarks/QC/2_deduplicated/{sample}.txt"
    resources:
        mem=config["simplejob_memory"],
        java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
    shell:
        """
        ( 
        if [ '{params.run_step}' == 't' ]; then
            clumpify.sh \
                {params.inputs} \
                {params.outputs} \
                overwrite=true \
                dedupe=t \
                dupesubs={params.dupesubs} \
                optical={params.only_optical} \
                threads={threads} \
                pigz=t unpigz=t \
                -Xmx{resources.java_mem}G
        else
            echo 'NOTE: Skipping step and hard linking files instead.'
            mkdir -p "{output}"
            for group in {params.pairs};
            do
                IFS=',' read -ra items <<< "$group"
                ln ${{items[0]}} ${{items[1]}}
            done
        fi
        ) 1>{log} 2>&1
        """


####
#### Filtered
####
PROCESSED_STEPS.append("3_quality_filtered")

rule apply_quality_filter:
    input:
        reads=lambda wc: get_output_fastq(wc, "2_deduplicated"),
        adapters=ancient(config["preprocess_adapters"]),
    output:
        reads=temp(lambda wc: get_output_fastq(wc, "3_quality_filtered")),
        stats="{sample}/logs/{sample}_quality_filtering_stats.txt",
    params:
        inputs =lambda wc, input:  io_params_for_tadpole(input.reads,  "in"),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, "out"),
        pairs=lambda wc, input, output: " ".join(",".join(x) for x in list(zip(input.reads, output.reads))),
        run_step=lambda wc: "t" if check_bool(wc, "Quality_filter") else "f",
        ref=(
            "ref=%s" % config["preprocess_adapters"]
            if (config["preprocess_adapters"] is not None)
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
            if PAIRED_END and config["error_correction_overlapping_pairs"]
            else "f"
        ),
        minlength=config["preprocess_minimum_passing_read_length"],
        minbasefrequency=config["preprocess_minimum_base_frequency"],
        # we require the user to reformat to R1 and R2, non-interleaved files
        interleaved="f",
        maxns=config["preprocess_max_ns"],
        prealloc=config["preallocate_ram"],
    log:
        "{sample}/logs/QC/3_quality_filtered.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    benchmark:
        "logs/benchmarks/QC/3_quality_filtered/{sample}.txt"
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
    shell:
        """
        ( 
        if [ '{params.run_step}' == 't' ]; then
            bbduk.sh {params.inputs} \
                {params.ref} \
                interleaved={params.interleaved} \
                {params.outputs} \
                stats={output.stats} \
                overwrite=true \
                qout=33 \
                trd=t \
                {params.hdist} \
                {params.k} \
                {params.ktrim} \
                {params.mink} \
                trimq={params.trimq} \
                qtrim={params.qtrim} \
                threads={threads} \
                minlength={params.minlength} \
                maxns={params.maxns} \
                minbasefrequency={params.minbasefrequency} \
                ecco={params.error_correction_pe} \
                prealloc={params.prealloc} \
                pigz=t unpigz=t \
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


####
#### Contaminant References
####
# if there are no references, decontamination will be skipped
if len(config.get("contaminant_references", {}).keys()) > 0:
    PROCESSED_STEPS.append("4_decontamination")

    rule build_decontamination_db:
        input:
            ancient(config["contaminant_references"].values()),
        output:
            "ref/genome/1/summary.txt",
        params:
            k=config["contaminant_kmer_length"],
            refs_in=" ".join(
                [
                    "ref_%s=%s" % (n, fa)
                    for n, fa in config["contaminant_references"].items()
                ]
            ),
        log:
            "logs/QC/4_decontamination_build_db.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["large_threads"]
        resources:
            mem=config["large_memory"],
            java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        shell:
            "bbsplit.sh"
            " -Xmx{resources.java_mem}G "
            " {params.refs_in} "
            " threads={threads}"
            " k={params.k}"
            " local=t "
            " &> {log}"

    rule run_decontamination:
        input:
            reads=lambda wc: get_output_fastq(wc, "3_quality_filtered"),
            db="ref/genome/1/summary.txt",
        output:
            reads=temp(lambda wc: get_output_fastq(wc, "4_decontamination")),
            stats="{sample}/sequence_quality_control/4_decontamination_reference_stats.txt",
            contaminant_folder=directory("{sample}/sequence_quality_control/4_contaminants")
        params:
            inputs =lambda wc, input:  io_params_for_tadpole(input.reads,  "in"),
            outputs=lambda wc, output: io_params_for_tadpole(output.reads, "outu"),
            pairs=lambda wc: " ".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step=lambda wc: "t" if check_bool(wc, "Remove_contaminants") else "f",
            paired="true" if PAIRED_END else "false",
            maxindel=config["contaminant_max_indel"],
            minratio=config["contaminant_min_ratio"],
            minhits=config["contaminant_minimum_hits"],
            ambiguous=config["contaminant_ambiguous"],
            k=config["contaminant_kmer_length"],
        log:
            "{sample}/logs/QC/4_decontamination.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["large_threads"]
        benchmark:
            "logs/benchmarks/QC/4_decontamination/{sample}.txt"
        resources:
            mem=config["large_memory"],
            java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        shell:
            """
            ( 
            if [ '{params.run_step}' == 't' ]; then
                bbsplit.sh \
                    {params.inputs} \
                    {params.outputs} \
                    basename={output.contaminant_folder}/%_R#.fastq.gz \
                    maxindel={params.maxindel} \
                    minratio={params.minratio} \
                    minhits={params.minhits} \
                    ambiguous={params.ambiguous} \
                    refstats={output.stats} \
                    threads={threads} \
                    k={params.k} \
                    local=t \
                    machineout=t \
                    pigz=t unpigz=t ziplevel=9 \
                    -Xmx{resources.java_mem}G
            else
                echo 'NOTE: Skipping step and hard linking files instead.'
                mkdir -p "{output.reads}"
                mkdir -p "{output.contaminant_folder}"
                touch "{output.stats}"
                for group in {params.pairs};
                do
                    IFS=',' read -ra items <<< "$group"
                    ln ${{items[0]}} ${{items[1]}}
                done
            fi
            ) 1>{log} 2>&1
            """

#                seal.sh \
#                    {params.inputs} \
#                    {params.outputs} \
#                    pattern={output.contaminant_folder}/%.fastq.gz \
#                    minhits={params.minhits} \
#                    ambiguous={params.ambiguous} \
#                    refstats={output.stats} \
#                    threads={threads} \
#                    k={params.k} \
#                    pigz=t unpigz=t ziplevel=9 \
#                    -Xmx{resources.java_mem}G


####
#### QC
####
PROCESSED_STEPS.append("5_QC")

localrules:
    qcreads,

rule qcreads:
    input:
        (
            lambda wc: get_output_fastq(wc, "4_decontamination")
            if len(config.get("contaminant_references", {}).keys()) > 0
            else lambda wc: get_output_fastq(wc, "3_quality_filtered")
        ),
    output:
        temp(lambda wc: get_output_fastq(wc, "5_QC")),
    params:
        inputs=lambda wc: get_ribosomal_rna_input(wc),
        outputs=lambda wc: get_output_fastq(wc, "5_QC"),
    run:
        import shutil, os
        import pandas as pd
        
        os.makedirs(str(output), exist_ok=True)
        for i in range(len(params.outputs)):
            with open(params.outputs[i], "wb") as outFile:
                with open(params.inputs['clean_reads'][i], "rb") as infile1:
                    shutil.copyfileobj(infile1, outFile)
                    if hasattr(params.inputs, "rrna_reads"):
                        with open(params.inputs['rrna_reads'][i], "rb") as infile2:
                            shutil.copyfileobj(infile2, outFile)


def get_qc_and_lr_reads(wildcards):
    """
    Returns a named dictionary of paths to QC'd reads. If R1, R2, and LR are all
    provided for a sample, the dictionary also contains the path to the long
    reads file for scaffolding. This can be unpacked by a rule.
    """
    inputs = {
        "qc_reads": get_output_fastq(wildcards, "5_QC")
    }

    sample_info = sampleTable.loc[wildcards.sample,].dropna()

    if (
        "Reads_raw_R1" in sample_info
        and "Reads_raw_R2" in sample_info
        and "Reads_raw_Long" in sample_info
    ):
        inputs["long_reads_for_scaffolding"] = sample_info["Reads_raw_Long"]

    return inputs


rule copy_qc_reads:
    input:
        unpack(get_qc_and_lr_reads),
    output:
        directory("QC/reads/{sample}"),
    params:
        sample=lambda wc: wc.sample,
    run:
        import shutil, os
        from glob import glob

        os.makedirs(str(output[0]), exist_ok=True)

        for f_in in input.qc_reads:
            shutil.copy(f_in, str(output[0]))

        if hasattr(input, "long_reads_for_scaffolding"):
            shutil.copy(input.long_reads_for_scaffolding, os.path.join(str(output[0]), params.sample + "_LR.fastq.gz"))


def get_quality_controlled_path(wildcards):
    """
    Gets path to quality controlled reads.
     - Just short hand for the final filtering step.
    """
    return("QC/reads/{sample}")

def get_quality_controlled_reads(wildcards, as_dict=False):
    if as_dict:
        files = {}
        for fraction in get_fractions(wildcards.sample):
            files[fraction].append(expand(
                "QC/reads/{sample}/{sample}_{fraction}.fastq.gz",
                sample=wildcards.sample,
                fraction=fraction,
            ))
    else:
        files = expand(
            "QC/reads/{sample}/{sample}_{fraction}.fastq.gz",
           sample=wildcards.sample,
            fraction=get_fractions(wildcards.sample),
        )
    return(files)


####
#### STATS
####

rule get_read_stats:
    input:
        "{sample}/sequence_quality_control/{step}",
    output:
        "{sample}/sequence_quality_control/read_stats/{step}.zip",
        read_counts=temp(
            "{sample}/sequence_quality_control/read_stats/{step}_read_counts.tsv"
        ),
    priority: 30
    log:
        "{sample}/logs/QC/read_stats/{step}.log",
    # conda:
    #     "../envs/required_packages.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        java_mem=int(config["simplejob_memory"] * JAVA_MEM_FRACTION),
    params:
        inputs=lambda wc: get_output_fastq(wc, wc.step),
        folder=lambda wc, output: os.path.splitext(output[0])[0],
        single_end_file=(
            "{sample}/sequence_quality_control/{sample}_{step}_se.fastq.gz"
        ),
    script:
        "../scripts/get_read_stats.py"


rule calculate_insert_size:
    input:
        "{sample}/sequence_quality_control/5_QC",
    output:
        ihist=(
            "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt"
        ),
        read_length=(
            "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        ),
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
    conda:
        "../envs/required_packages.yaml"
    log:
        "{sample}/logs/QC/stats/calculate_insert_size.log",
    params:
        kmer=config["merging_k"],
        extend2=config["merging_extend2"],
        flags="loose ecct",
        minprob=config.get("bbmerge_minprob", "0.8"),
        inputs=lambda wc: io_params_for_tadpole(get_quality_controlled_reads(wc)),
        check_insert=lambda wc: "t" if check_paired(wc.sample) else "f",
    shell:
        """
        ( 
        readlength.sh {params.inputs} out={output.read_length}
        
        if [ '{params.check_insert}' == 't' ]; then
            bbmerge.sh \
                -Xmx{resources.java_mem}G \
                threads={threads} \
                {params.inputs} \
                {params.flags} k={params.kmer} \
                extend2={params.extend2} \
                ihist={output.ihist} merge=f \
                mininsert0=35 minoverlap0=8 \
                prealloc=t prefilter=t \
                minprob={params.minprob}
        else
            echo 'NOTE: Skipping Inster Size Hist step. Creating empty output file.'
            touch {output.ihist}
        fi
        ) 1>{log} 2>&1
        """


localrules:
    combine_read_length_stats,
    combine_insert_stats,

rule combine_read_length_stats:
    input:
        expand(
            "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
            sample=SAMPLES,
        ),
    output:
        "stats/read_length_stats.tsv",
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()

        for length_file in input:
            sample = length_file.split(os.path.sep)[0]
            data = parse_comments(length_file)
            data = pd.Series(data)[
                ["Reads", "Bases", "Max", "Min", "Avg", "Median", "Mode", "Std_Dev"]
            ]
            stats[sample] = data

        stats.to_csv(output[0], sep="\t")


rule combine_insert_stats:
    input:
        expand(
            "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",
            sample=SAMPLES,
        ),
    output:
        "stats/insert_stats.tsv",
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()

        for insert_file in input:
            sample = insert_file.split(os.path.sep)[0]
            if os.path.isfile(insert_file) and os.path.getsize(insert_file) == 0:
                continue
            data = parse_comments(insert_file)
            data = pd.Series(data)[
                ["Mean", "Median", "Mode", "STDev", "PercentOfPairs"]
            ]
            stats[sample] = data

        stats.T.to_csv(output[0], sep="\t")


localrules:
    combine_read_counts,
    write_read_counts,

rule write_read_counts:
    input:
        read_count_files=expand(
            "{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv",
            step=PROCESSED_STEPS,
        ),
    output:
        read_stats="{sample}/sequence_quality_control/read_stats/read_counts.tsv",
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
            "{sample}/sequence_quality_control/read_stats/read_counts.tsv",
            sample=SAMPLES,
        ),
    output:
        "stats/read_counts.tsv",
    run:
        from utils.io import pandas_concat

        pandas_concat(list(input), output[0], sep="\t", index_col=[0, 1], axis=0)


rule build_qc_report:
    input:
        zipfiles_QC=expand(
            "{sample}/sequence_quality_control/read_stats/5_QC.zip", sample=SAMPLES
        ),
        read_counts="stats/read_counts.tsv",
        read_length_stats=(
            ["stats/read_length_stats.tsv", "stats/insert_stats.tsv"]
            if PAIRED_END
            else "stats/read_length_stats.tsv"
        ),
    output:
        report="reports/QC_report.html",
    log:
        "logs/QC/report.log",
    params:
        min_quality=config["preprocess_minimum_base_quality"],
        samples=SAMPLES,
    conda:
        "../envs/report.yaml"
    script:
        "../report/qc_report.py"


####
#### Done
####

rule finalize_sample_qc:
    input:
        reads="QC/reads/{sample}",
        reads_stats_zip=expand(
            "{{sample}}/sequence_quality_control/read_stats/{step}.zip",
            step=PROCESSED_STEPS,
        ),
        read_length_hist=(
            "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        ),
        report="reports/QC_report.html",
    output:
        flag=touch("{sample}/sequence_quality_control/finished_QC"),
