import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings
from copy import deepcopy


def get_preprocessing_steps(config):
    preprocessing_steps = ["QC"]
    if config.get("error_correction_before_assembly", True):
        preprocessing_steps.append("errorcorr")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)


####
#### Normalize Reads
####
def normalize_reads_command(inputs, outputs, outdir, pairs, histin, histout, run_step, k, target, mindepth, threads, resources):
    if run_step == 't':
        cmd = f"""
        bbnorm.sh \\
            {inputs} \\
            {outputs} \\
            tmpdir={resources.tmpdir} \\
            tossbadreads=t \\
            hist={histin} \\
            histout={histout} \\
            mindepth={mindepth} \\
            k={k} \\
            target={target} \\
            prefilter=t \\
            threads={threads} \\
            -Xmx{resources.java_mem}G
        """
    else:
        cmd = f"""
        echo 'Skipping step, hard linking files instead.'
        mkdir -p "{outdir}"
        touch "{histin}"
        touch "{histout}"
        IFS=';' read -ra groups <<< '{pairs}'
        for group in "${{groups[@]}}";
        do
            IFS=',' read -ra items <<< "$group"
            cp "${{items[0]}}" "${{items[1]}}"
        done
        """
    return(cmd)


rule normalize_reads_PE:
    input:
        reads=[
            "{sample}/sequence_quality_control/{sample}_R1.fastq.gz",
            "{sample}/sequence_quality_control/{sample}_R2.fastq.gz"
        ],
    output:
        reads=temp([
            "{sample}/assembly/reads/1_normalize_reads_R1.fastq.gz",
            "{sample}/assembly/reads/1_normalize_reads_R2.fastq.gz"
        ]),
        histin ="{sample}/assembly/reads/1_normalize_reads_PE.histogram_before_normalization.tsv.gz",
        histout="{sample}/assembly/reads/1_normalize_reads_PE.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wc, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wc, "Normalize_reads_before_assembly") else "f",
            k=config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
            target=config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/1_normalize_reads.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/1_normalize_reads/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule normalize_reads_SE:
    input:
        reads=[
            "{sample}/sequence_quality_control/{sample}_SE.fastq.gz",
        ],
    output:
        reads=temp([
            "{sample}/assembly/reads/1_normalize_reads_SE.fastq.gz",
        ]),
        histin ="{sample}/assembly/reads/1_normalize_reads_SE.histogram_before_normalization.tsv.gz",
        histout="{sample}/assembly/reads/1_normalize_reads_SE.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wc, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wc, "Normalize_reads_before_assembly") else "f",
            k=config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
            target=config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/1_normalize_reads.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/1_normalize_reads/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule normalize_reads_LR:
    input:
        reads=[
            "{sample}/sequence_quality_control/{sample}_LR.fastq.gz",
        ],
    output:
        reads=temp([
            "{sample}/assembly/reads/1_normalize_reads_LR.fastq.gz",
        ]),
        histin ="{sample}/assembly/reads/1_normalize_reads_LR.histogram_before_normalization.tsv.gz",
        histout="{sample}/assembly/reads/1_normalize_reads_LR.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wc, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wc, "Normalize_reads_before_assembly") else "f",
            k=config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
            target=config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/1_normalize_reads.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/1_normalize_reads/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """



####
#### Error Correction
####
def error_correction_command(inputs, outputs, outdir, pairs, run_step, 
                            prefilter, minprob, tossdepth, tossjunk, lowdepthfraction, 
                            aggressive, shave, threads, resources):
    if run_step == 't':
        cmd = f"""
            tadpole.sh \\
                prefilter={prefilter} \\
                prealloc=1 \\
                {inputs} \\
                {outputs} \\
                mode=correct \\
                aggressive={aggressive} \\
                tossjunk={tossjunk} \\
                lowdepthfraction={lowdepthfraction} \\
                tossdepth={tossdepth} \\
                merge=t \\
                shave={shave} \\
                rinse={shave} \\
                threads={threads} \\
                pigz=t \\
                unpigz=t \\
                ecc=t \\
                ecco=t \\
                -Xmx{resources.java_mem}G
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

rule error_correction_PE:
    input:
        reads=rules.normalize_reads_PE.output.reads,
    output:
        reads=temp([
            "{sample}/assembly/reads/2_error_correction_R1.fastq.gz",
            "{sample}/assembly/reads/2_error_correction_R2.fastq.gz"
        ]),
    params:
        command = lambda wc, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Error_correction_before_assembly") else "f",
            prefilter=2,  # Ignore kmers with less than 2 occurance
            minprob=config["error_correction_minprob"],
            tossdepth=config["error_correction_minimum_kmer_depth"],
            tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
            lowdepthfraction=config["error_correction_lowdepth_fraction"],
            aggressive=config["error_correction_aggressive"],
            shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/2_error_correction.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/2_error_correction/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule error_correction_SE:
    input:
        reads=rules.normalize_reads_SE.output.reads,
    output:
        reads=temp([
            "{sample}/assembly/reads/2_error_correction_SE.fastq.gz",
        ]),
    params:
        command = lambda wc, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Error_correction_before_assembly") else "f",
            prefilter=2,  # Ignore kmers with less than 2 occurance
            minprob=config["error_correction_minprob"],
            tossdepth=config["error_correction_minimum_kmer_depth"],
            tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
            lowdepthfraction=config["error_correction_lowdepth_fraction"],
            aggressive=config["error_correction_aggressive"],
            shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/2_error_correction.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/2_error_correction/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule error_correction_LR:
    input:
        reads=rules.normalize_reads_LR.output.reads,
    output:
        reads=temp([
            "{sample}/assembly/reads/2_error_correction_LR.fastq.gz",
        ]),
    params:
        command = lambda wc, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"{wc.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wc, "Error_correction_before_assembly") else "f",
            prefilter=2,  # Ignore kmers with less than 2 occurance
            minprob=config["error_correction_minprob"],
            tossdepth=config["error_correction_minimum_kmer_depth"],
            tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
            lowdepthfraction=config["error_correction_lowdepth_fraction"],
            aggressive=config["error_correction_aggressive"],
            shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
            threads=threads,
            resources=resources
        )
    log:
        "{sample}/logs/assembly/pre_process/2_error_correction.log",
    benchmark:
        "{sample}/benchmarks/assembly/pre_process/2_error_correction/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],    
        java_mem=int(config["large_memory"] * JAVA_MEM_FRACTION),
        time=config["large_runtime"], 
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


def get_pre_processed_reads(wildcards, as_dict=False):
    if as_dict:
        files = {}
        for fraction in get_fractions(wildcards.sample):
            files[fraction] = f"{wildcards.sample}/assembly/reads/2_error_correction_{fraction}.fastq.gz"
    else:
        files = expand(
            "{sample}/assembly/reads/2_error_correction_{fraction}.fastq.gz",
            sample=wildcards.sample,
            fraction=get_fractions(wildcards.sample),
        )
    return(files)



####
#### Assembly
####

def assembly_command(wildcards, input, output, threads, resources):
    """
    Builds the assembly command depending on the type of data we have.
    """
    assembler = sampleTable.loc[wildcards.sample, 'Assembler']
    output_dir = f"{wildcards.sample}/assembly/assembly"
    
    # SPADES (Short Reads + long reads for scaffolding)
    if assembler.startswith('spades'):
        s_num = 1 # Keep track of how many SE files we have
        
        reads = ""
        # Check for named inputs 'r1'/'r2' or 'se'
        if hasattr(input, "R1"): 
            reads += f" -1 {input.R1} -2 {input.R2} "
        elif hasattr(input, "SE"):
            reads += f" --s{s_num} {input.SE} "
            s_num+=1
        else:
            raise ValueError(f"No trimmed short reads found for assembler '{assembler}' and sample '{wildcards.sample}'.")
        
        # Add long reads for SPADES hybrid assembly
        if hasattr(input, "LR"):
            # HQ PacBio and Nanopore should be set as single end reads.
            # See: https://ablab.github.io/spades/running.html
            m = {"spades-pacbio-raw":   "--pacbio",
                 "spades-pacbio-corr":  "--pacbio",
                 "spades-pacbio-hq":    f"--s{s_num}",
                 "spades-nanopore-raw": "--nanopore",
                 "spades-nanopore-corr":"--nanopore",
                 "spades-nanopore-hq":  f"--s{s_num}"}
            reads += f" {m[assembler]} {input.LR} "
        
        k = config.get("spades_k", SPADES_K)
        extra = config.get("spades_extra", '')
        sequences="scaffolds" if config["spades_use_scaffolds"] else "contigs"
        
        # If we see spades has already run, we can continue to save time.
        if not os.path.exists(f"{output_dir}/params.txt"):
            cmd = f"""
            rm -fr "{output_dir}"
            
            spades.py \\
                {reads} \\
                -o {output_dir} \\
                --meta \\
                --only-assembler \\
                -k {k} \\
                --checkpoints last \\
                --threads {threads} \\
                --memory {resources.mem} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{sequences}.fasta" > {output}
            """
        else:
            cmd = f"""
            spades.py \\
                -o {output_dir} \\
                --restart-from last \\
                -k {k} \\
                --threads {threads} \\
                --memory {resources.mem} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{sequences}.fasta" > {output}
            """
    
    
    # MEGAHIT (Short Reads)
    elif assembler.startswith('megahit'):
        reads = ""
        # Check for named inputs 'r1'/'r2' or 'se'
        if hasattr(input, "R1"):
            reads += f"-1 {input.R1} -2 {input.R2}"
        elif hasattr(input, "SE"):
            reads += f"-r {input.SE}"
        else:
            raise ValueError(f"No trimmed short reads found for assembler '{assembler}' and sample '{sample_id}'.")
        
        min_count=config.get("megahit_min_count", MEGAHIT_MIN_COUNT),
        k_min=config.get("megahit_k_min", MEGAHIT_K_MIN),
        k_max=config.get("megahit_k_max", MEGAHIT_K_MAX),
        k_step=config.get("megahit_k_step", MEGAHIT_K_STEP),
        merge_level=config.get("megahit_merge_level", MEGAHIT_MERGE_LEVEL),
        prune_level=config.get("megahit_prune_level", MEGAHIT_PRUNE_LEVEL),
        low_local_ratio=config["megahit_low_local_ratio"],
        min_contig_len=config["minimum_contig_length"],
        assembly_params = {
            "default": "",
            "meta-sensitive": "--presets meta-sensitive",
            "meta-large": " --presets meta-large",
        }
        preset=assembly_params[config["megahit_preset"]],
        extra = config.get("megahit_extra", '')
        
        # If we see megahit has already run, we can continue to save time.
        if not os.path.exists(f"{output_dir}/options.json"):
            cmd = f"""
            rm -fr "{output_dir}"
            
            megahit \\
                {reads} \\
                --out-dir {output_dir} \\
                --out-prefix {wildcards.sample}_prefilter \\
                --tmp-dir {resources.tmpdir} \\
                --num-cpu-threads {threads} \\
                --k-min {k_min[0]} \\
                --k-max {k_max[0]} \\
                --k-step {k_step[0]} \\
                --min-contig-len {min_contig_len[0]} \\
                --min-count {min_count[0]} \\
                --merge-level {merge_level[0]} \\
                --prune-level {prune_level[0]} \\
                --low-local-ratio {low_local_ratio[0]} \\
                --memory {resources.mem}000000 \\
                {preset[0]} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{wildcards.sample}_prefilter.contigs.fa" > {output}
            """
        else:
            cmd = f"""
            megahit \\
                --out-dir {output_dir} \\
                --continue
            
            seqkit sort -l -r -w 0 "{output_dir}/{wildcards.sample}_prefilter.contigs.fa" > {output}
            """
    
    
    # FLYE (Long Reads)
    elif assembler.startswith('flye'):
        extra = config.get("flye_extra", '')
        
        m = {"flye-pacbio-raw":   "--pacbio-raw", 
             "flye-pacbio-corr":  "--pacbio-corr",
             "flye-pacbio-hq":    "--pacbio-hq",
             "flye-nanopore-raw": "--nano-raw",
             "flye-nanopore-corr":"--nano-corr",
             "flye-nanopore-hq":  "--nano-hq"}
        reads = f"{m[assembler]} {input.LR}"

        # If we see flye has already run, we can continue to save time.
        if not os.path.exists(f"{output_dir}/params.json"):
            cmd = f"""
            rm -fr "{output_dir}"
            
            flye \\
                {reads} \\
                --out-dir {output_dir} \\
                --meta \\
                --threads {threads} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/assembly.fasta" > {output}
            """
        else:
            cmd = f"""
            flye \\
                {reads} \\
                --out-dir {output_dir} \\
                --meta \\
                --threads {threads} {extra} \\
                --resume
            
            seqkit sort -l -r -w 0 "{output_dir}/assembly.fasta" > {output}
            """
    
    
    # metaMDBG (Long Reads)
    elif assembler.startswith('metamdbg'):
        extra = config.get("metamdbg_extra", '')
        
        m = {"metamdbg-pacbio-hq":  "--in-hifi",
             "metamdbg-nanopore-hq":"--in-ont"}
        reads = f"{m[assembler]} {input.LR}"
        
        # If metaMDBG has already run, it should resume automatically.
        cmd = f"""
        metaMDBG asm \\
            {reads} \\
            --out-dir {output_dir} \\
            --skip-correction \\
            --threads {threads} {extra}
        
        zcat "{output_dir}/contigs.fasta.gz" | sed -e 's/ .*circular=/_circular_/' | seqkit sort -l -r -w 0 > {output}
        """
    
    
    # Dont recognize assembler
    else:
        raise ValueError(f"Unknown assembler '{assembler}'.")
    
    return(cmd)


rule run_assembly:
    input:
        unpack(lambda wc: get_pre_processed_reads(wc, as_dict=True)),
    output:
        "{sample}/assembly/assembly/{sample}_raw_contigs.fasta"
    params:
        command = lambda wildcards, input, output, threads, resources: assembly_command(
            wildcards, input, output, threads, resources
        ),
    log:
        "{sample}/logs/assembly.log",
    benchmark:
        "{sample}/benchmarks/assembly/{sample}.txt"
    conda:
        "../envs/assembly.yaml"
    threads: config["assembly_threads"]
    resources:
        mem=config["assembly_memory"],
        time_min=60 * config["assembly_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule rename_contigs:
    input:
        "{sample}/assembly/assembly/{sample}_raw_contigs.fasta",
    output:
        fasta="{sample}/assembly/assembly/{sample}_prefilter_contigs.fasta",
        mapping_table="{sample}/assembly/assembly/old2new_contig_names.tsv",
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


def align_reads_command(wildcards, input, output, threads, resources):
    """
    Map any combination of reads against a reference using minimap2.
    """
    
    cmd = ""
    
    ## Map short reads
    cmd_sr = ""
    if hasattr(input, "R1"):
        cmd_sr = f"minimap2 -t {threads} -ax sr {input.target} {input.R1} {input.R2}"
    elif hasattr(input, "SE"):
        cmd_sr = f"minimap2 -t {threads} -ax sr {input.target} {input.SE}"
    
    ## Map long reads
    # pacbio-raw:      PacBio regular CLR reads (<20% error)
    # pacbio-corr:     PacBio reads that were corrected with other methods (<3% error)
    # pacbio-hq:       PacBio HiFi reads (<1% error)
    # nanopore-raw:    ONT regular reads, pre-Guppy5 (<20% error)
    # nanopore-corr:   ONT reads that were corrected with other methods (<3% error)
    # nanopore-hq:     ONT high-quality reads (<1% error)
    cmd_lr = ""
    if hasattr(input, "LR"):
        assembler = sampleTable.loc[wildcards.sample, 'Assembler']
        if assembler.endswith("-pacbio-raw") | assembler.endswith("-pacbio-corr"):
            cmd_lr = f"minimap2 -t {threads} -ax map-pb {input.target} {input.LR}"
        
        elif assembler.endswith("-nanopore-raw") | assembler.endswith("-nanopore-corr"):
            cmd_lr = f"minimap2 -t {threads} -ax map-ont {input.target} {input.LR}"
        
        elif assembler.endswith("-pacbio-hq"):
            cmd_lr = f"minimap2 -t {threads} -ax map-hifi {input.target} {input.LR}"
        
        elif assembler.endswith("-nanopore-hq"):
            # See: https://github.com/lh3/minimap2/issues/1127
            cmd_lr = f"minimap2 -t {threads} -ax lr:hq {input.target} {input.LR}"
        
        # Unknown error rate
        else:
            cmd_lr = f"minimap2 -t {threads} -ax map-ont {input.target} {input.LR}"
   
    # Check if we have SR+LR (need to map separatly and merge) or SR OR LR
    if cmd_sr and cmd_lr:
        cmd = f"({cmd_sr} && {cmd_lr} | grep -v '^@') | samtools sort > {output}"
    elif cmd_sr and not cmd_lr:
        cmd = f"{cmd_sr} | samtools sort > {output}"
    else:
        cmd = f"{cmd_lr} | samtools sort > {output}"
    
    return(cmd)


if config["filter_contigs"]:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            unpack(lambda wc: get_pre_processed_reads(wc, as_dict=True)),
            target=rules.rename_contigs.output.fasta,
        output:
            bam=temp("{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam"),
        params:
            command = lambda wildcards, input, output, threads, resources: align_reads_command(
                wildcards, input, output, threads, resources
            ),
        benchmark:
            "{sample}/benchmarks/assembly/post_process/align_reads_to_prefiltered_contigs.txt",
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log",
        conda:
            "../envs/minimap.yaml"
        threads: config["simplejob_threads"]
        resources:
            mem=config["simplejob_memory"],
            time=config["simplejob_runtime"],
        shell:
            """
            ({params.command}) >{log} 2>&1
            """


    rule pileup_prefilter:
        input:
            fasta="{sample}/assembly/assembly/{sample}_prefilter_contigs.fasta",
            bam="{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam",
        output:
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        params:
            pileup_secondary="t",
            minmapq=config["minimum_map_quality"],
        benchmark:
            "{sample}/benchmarks/assembly/post_process/pilup_prefilter_contigs.log",
        log:
            "{sample}/logs/assembly/post_process/pilup_prefilter_contigs.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["simplejob_threads"]
        resources:
            mem=config["simplejob_memory"],
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
            fasta="{sample}/assembly/assembly/{sample}_prefilter_contigs.fasta",
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        output:
            fasta="{sample}/assembly/assembly/{sample}_final_contigs.fasta",
            removed_names="{sample}/assembly/assembly/{sample}_discarded_contigs.fasta",
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
            "{sample}/assembly/assembly/{sample}_prefilter_contigs.fasta",
        output:
            "{sample}/assembly/assembly/{sample}_final_contigs.fasta",
        shell:
            "cp {input} {output}"


localrules:
    finalize_contigs,


rule finalize_contigs:
    input:
        "{sample}/assembly/assembly/{sample}_final_contigs.fasta",
    output:
        "{sample}/assembly/{sample}.fasta",
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
        unpack(lambda wc: get_pre_processed_reads(wc, as_dict=True)),
        target="{sample}/assembly/{sample_contigs}.fasta",
    output:
        bam=temp("{sample_contigs}/sequence_alignment/{sample}.bam"),
    params:
        command = lambda wildcards, input, output, threads, resources: align_reads_command(
            wildcards, input, output, threads, resources
        ),
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/align_reads_to_filtered_contigs/{sample}_to_{sample_contigs}.txt",
    log:
        "{sample_contigs}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log",
    conda:
        "../envs/minimap.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


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
