import os
import re
import sys
import warnings
from glob import glob
from copy import deepcopy



####
#### Normalize Reads
####
def normalize_reads_command(inputs, outputs, outdir, pairs, histin, histout, run_step, k, target, mindepth, threads, resources):
    if run_step == 't':
        cmd = f"""
        bbnorm.sh \\
            {inputs} \\
            {outputs} \\
            tmpdir={outdir}/tmp \\
            tossbadreads=t \\
            hist={histin} \\
            histout={histout} \\
            mindepth={mindepth} \\
            k={k} \\
            target={target} \\
            prefilter=t \\
            threads={threads} \\
            -Xmx{resources.java_mem}M
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
            "samples/{sample}/sequence_quality_control/{sample}_R1.fastq.gz",
            "samples/{sample}/sequence_quality_control/{sample}_R2.fastq.gz"
        ],
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/1_normalize_reads_R1.fastq.gz",
            "samples/{sample}/assembly/reads/1_normalize_reads_R2.fastq.gz"
        ]),
        histin ="samples/{sample}/assembly/reads/1_normalize_reads_PE.histogram_before_normalization.tsv.gz",
        histout="samples/{sample}/assembly/reads/1_normalize_reads_PE.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wildcards, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wildcards, "Normalize_reads_before_assembly") else "f",
            k=config["normalization_kmer_length"],
            target=config["normalization_target_depth"],
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        ),
        tmp="samples/{sample}/assembly/reads/tmp",
    log:
        "logs/samples/{sample}/assembly/reads/1_normalize_reads_PE.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/1_normalize_reads_PE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "normalize_reads", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule normalize_reads_SE:
    input:
        reads=[
            "samples/{sample}/sequence_quality_control/{sample}_SE.fastq.gz",
        ],
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/1_normalize_reads_SE.fastq.gz",
        ]),
        histin ="samples/{sample}/assembly/reads/1_normalize_reads_SE.histogram_before_normalization.tsv.gz",
        histout="samples/{sample}/assembly/reads/1_normalize_reads_SE.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wildcards, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wildcards, "Normalize_reads_before_assembly") else "f",
            k=config["normalization_kmer_length"],
            target=config["normalization_target_depth"],
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/assembly/reads/1_normalize_reads_SE.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/1_normalize_reads_SE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "normalize_reads", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule normalize_reads_LR:
    input:
        reads=[
            "samples/{sample}/sequence_quality_control/{sample}_LR.fastq.gz",
        ],
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/1_normalize_reads_LR.fastq.gz",
        ]),
        histin ="samples/{sample}/assembly/reads/1_normalize_reads_LR.histogram_before_normalization.tsv.gz",
        histout="samples/{sample}/assembly/reads/1_normalize_reads_LR.histogram_after_normalization.tsv.gz",
    params:
        command = lambda wildcards, input, output, threads, resources: normalize_reads_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            histin=output.histin,
            histout=output.histout,
            run_step="t" if check_bool(wildcards, "Normalize_reads_before_assembly") else "f",
            k=config["normalization_kmer_length"],
            target=config["normalization_target_depth"],
            mindepth=config["normalization_minimum_kmer_depth"],
            threads=threads,
            resources=resources
        )
    log:
        "logs/samples/{sample}/assembly/reads/1_normalize_reads_LR.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/1_normalize_reads_LR.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "normalize_reads", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "normalize_reads", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
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

rule error_correction_PE:
    input:
        reads=rules.normalize_reads_PE.output.reads,
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/2_error_correction_R1.fastq.gz",
            "samples/{sample}/assembly/reads/2_error_correction_R2.fastq.gz"
        ]),
    params:
        command = lambda wildcards, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wildcards, "Error_correction_before_assembly") else "f",
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
        "logs/samples/{sample}/assembly/reads/2_error_correction_PE.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/2_error_correction_PE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "error_correction", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule error_correction_SE:
    input:
        reads=rules.normalize_reads_SE.output.reads,
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/2_error_correction_SE.fastq.gz",
        ]),
    params:
        command = lambda wildcards, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wildcards, "Error_correction_before_assembly") else "f",
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
        "logs/samples/{sample}/assembly/reads/2_error_correction_SE.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/2_error_correction_SE.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "error_correction", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule error_correction_LR:
    input:
        reads=rules.normalize_reads_LR.output.reads,
    output:
        reads=temp([
            "samples/{sample}/assembly/reads/2_error_correction_LR.fastq.gz",
        ]),
    params:
        command = lambda wildcards, input, output, threads, resources: error_correction_command(
            inputs=io_params_for_tadpole(input.reads),
            outputs=io_params_for_tadpole(output.reads, key="out"),
            outdir=f"samples/{wildcards.sample}/assembly/reads",
            pairs=";".join(",".join(x) for x in list(zip(input.reads, output.reads))),
            run_step="t" if check_bool(wildcards, "Error_correction_before_assembly") else "f",
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
        "logs/samples/{sample}/assembly/reads/2_error_correction_LR.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/reads/2_error_correction_LR.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "error_correction", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "error_correction", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


def get_pre_processed_reads(wildcards, as_dict=False):
    if as_dict:
        files = {}
        for fraction in get_fractions(wildcards.sample):
                files[fraction] = f"samples/{wildcards.sample}/assembly/reads/2_error_correction_{fraction}.fastq.gz"
    else:
        files = expand(
                "samples/{sample}/assembly/reads/2_error_correction_{fraction}.fastq.gz",
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
    output_dir = f"{output.outdir}"

    #mem_gb = 1
    # To prevent lazy evaluation issues. Will be undefined or str during DAG construction, can only be calculated during rule execution.
    #if isinstance(resources.mem, (int, float)):
    #  mem_gb = resources.mem // 1024 # MB -> GB
    
    # SPADES (Short Reads + long reads for scaffolding)
    if assembler.startswith('spades'):
        s_num = 1 # Keep track of how many SE files we have
        
        reads = ""
        # Check for named inputs 'r1'/'r2' or 'se'
        if hasattr(input, "R1"): 
            reads += f"-1 {input.R1} -2 {input.R2} "
        elif hasattr(input, "SE"):
            reads += f"--s{s_num} {input.SE} "
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
        
        k = config["spades_k"]
        extra = config["spades_extra"]
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
                --memory ${{MEM_GB}} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{sequences}.fasta" > {output.fasta}
            """
        else:
            cmd = f"""
            spades.py \\
                -o {output_dir} \\
                --restart-from last \\
                -k {k} \\
                --threads {threads} \\
                --memory ${{MEM_GB}} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{sequences}.fasta" > {output.fasta}
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
        
        min_count=config["megahit_min_count"]
        k_min=config["megahit_k_min"]
        k_max=config["megahit_k_max"]
        k_step=config["megahit_k_step"]
        merge_level=config["megahit_merge_level"]
        prune_level=config["megahit_prune_level"]
        low_local_ratio=config["megahit_low_local_ratio"]
        min_contig_len=config["minimum_contig_length"]
        assembly_params = {
            "default": "",
            "meta-sensitive": "--presets meta-sensitive",
            "meta-large": " --presets meta-large",
        }
        preset=assembly_params[config["megahit_preset"]]
        extra = config["megahit_extra"]
        
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
                --k-min {k_min} \\
                --k-max {k_max} \\
                --k-step {k_step} \\
                --min-contig-len {min_contig_len} \\
                --min-count {min_count} \\
                --merge-level {merge_level} \\
                --prune-level {prune_level} \\
                --low-local-ratio {low_local_ratio} \\
                --memory ${{MEM_GB}}000000000 {preset} {extra}
            
            seqkit sort -l -r -w 0 "{output_dir}/{wildcards.sample}_prefilter.contigs.fa" > {output.fasta}
            """
        else:
            cmd = f"""
            megahit \\
                --out-dir {output_dir} \\
                --num-cpu-threads {threads} \\
                --memory ${{MEM_GB}}000000000 \\
                --continue
            
            seqkit sort -l -r -w 0 "{output_dir}/{wildcards.sample}_prefilter.contigs.fa" > {output.fasta}
            """
    
    
    # FLYE (Long Reads)
    elif assembler.startswith('flye'):
        extra = config["flye_extra"]
        
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
            
            seqkit sort -l -r -w 0 "{output_dir}/assembly.fasta" > {output.fasta}
            """
        else:
            cmd = f"""
            flye \\
                {reads} \\
                --out-dir {output_dir} \\
                --meta \\
                --threads {threads} {extra} \\
                --resume
            
            seqkit sort -l -r -w 0 "{output_dir}/assembly.fasta" > {output.fasta}
            """
    
    
    # metaMDBG (Long Reads)
    elif assembler.startswith('metamdbg'):
        extra = config["metamdbg_extra"]
        
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
        
        zcat "{output_dir}/contigs.fasta.gz" | sed -e 's/ .*circular=/_circular_/' | seqkit sort -l -r -w 0 > {output.fasta}
        """
    
    
    # Dont recognize assembler
    else:
        raise ValueError(f"Unknown assembler '{assembler}'.")
    
    return(cmd)


rule run_assembly:
    input:
        unpack(lambda wildcards: get_pre_processed_reads(wildcards, as_dict=True)),
    output:
        fasta="samples/{sample}/assembly/{sample}_raw_contigs.fasta",
	outdir=temp(directory("samples/{sample}/assembly/assembly")),
    params:
        command = lambda wildcards, input, output, threads, resources: assembly_command(
            wildcards, input, output, threads, resources
        ),
    log:
        "logs/samples/{sample}/assembly/assembly.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/assembly.tsv",
    conda:
        "../envs/assembly.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "run_assembly", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "run_assembly", "mem_mb"),
        mem_gb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "run_assembly", "mem_gb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "run_assembly", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "run_assembly", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "run_assembly", "account"),
    shell:
        """
        (MEM_GB={resources.mem_gb}; {params.command}) 1>{log} 2>&1
        """


rule rename_contigs:
    input:
        "samples/{sample}/assembly/{sample}_raw_contigs.fasta",
    output:
        fasta="samples/{sample}/assembly/{sample}_prefilter_contigs.fasta",
        mapping_table="samples/{sample}/assembly/{sample}_prefilter_old2new_contig_names.tsv",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "rename_contigs", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "rename_contigs", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "rename_contigs", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "rename_contigs", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "rename_contigs", "account"),
    log:
        "logs/samples/{sample}/assembly/rename_and_filter_size.log",
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
    
    cmd = "rm -fr {output}*; "
    
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
            ax = f"map-pb"
        elif assembler.endswith("-nanopore-raw") | assembler.endswith("-nanopore-corr"):
            ax = f"map-ont"
        elif assembler.endswith("-pacbio-hq"):
            ax = f"map-hifi"
        elif assembler.endswith("-nanopore-hq"):
            ax = f"lr:hq" # See: https://github.com/lh3/minimap2/issues/1127
        else:
            ax = f"map-ont"
        cmd_lr = f"minimap2 -t {threads} -ax {ax} {input.target} {input.LR}"
   
    # Check if we have SR+LR (need to map separatly and merge) or SR OR LR
    if cmd_sr and cmd_lr:
        cmd = f"{cmd} ({cmd_sr} && {cmd_lr} | grep -v '^@') | samtools sort"
    elif cmd_sr and not cmd_lr:
        cmd = f"{cmd} {cmd_sr} | samtools sort"
    else:
        cmd = f"{cmd} {cmd_lr} | samtools sort"
    
    # Set samtools sort temp file location
    cmd = f"{cmd} -T {output}.temp"
    
    # If output file provided, else will be printed to stdout
    if not output is None:
        cmd = f"{cmd} > {output}"
    
    return(cmd)


rule align_reads_to_prefilter_contigs:
    input:
        unpack(lambda wildcards: get_quality_controlled_reads(wildcards, as_dict=True)),
        target=rules.rename_contigs.output.fasta,
    output:
        bam=temp("samples/{sample}/assembly/{sample}_prefilter_contigs.bam"),
    params:
        command = lambda wildcards, input, output, threads, resources: align_reads_command(
            wildcards, input, output, threads, resources,
        ),
    benchmark:
        "benchmarks/samples/{sample}/assembly/align_reads_to_prefiltered_contigs.tsv",
    log:
        "logs/samples/{sample}/assembly/align_reads_to_prefiltered_contigs.log",
    conda:
        "../envs/minimap.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "mapping", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "account"),
    shell:
        """
        ({params.command}) 1>{log} 2>&1
        """


rule pileup_prefilter:
    input:
        fasta="samples/{sample}/assembly/{sample}_prefilter_contigs.fasta",
        bam="samples/{sample}/assembly/{sample}_prefilter_contigs.bam",
    output:
        covstats="samples/{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
    params:
        pileup_secondary="t",
        minmapq=config["minimum_map_quality"],
    benchmark:
        "benchmarks/samples/{sample}/assembly/pilup_prefilter_contigs.tsv",
    log:
        "logs/samples/{sample}/assembly/pilup_prefilter_contigs.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "pileup", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "account"),
    shell:
        """
        (
        pileup.sh \\
            ref={input.fasta} \\
            in={input.bam} \\
            threads={threads} \\
            -Xmx{resources.java_mem}M \\
            covstats={output.covstats} \\
            concise=t \\
            minmapq={params.minmapq} \\
            secondary={params.pileup_secondary}
        ) 1>{log} 2>&1
        """

rule filter_by_coverage:
    input:
        fasta="samples/{sample}/assembly/{sample}_prefilter_contigs.fasta",
        covstats="samples/{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
    output:
        fasta="samples/{sample}/assembly/{sample}_final_contigs.fasta",
        removed_names="samples/{sample}/assembly/{sample}_discarded_contigs.fasta",
    params:
        minc=config["minimum_average_coverage"],
        minp=config["minimum_percent_covered_bases"],
        minr=config["minimum_mapped_reads"],
        minl=config["minimum_contig_length"],
        trim=config["contig_trim_bp"],
    benchmark:
        "benchmarks/samples/{sample}/assembly/filter_by_coverage.tsv",
    log:
        "logs/samples/{sample}/assembly/filter_by_coverage.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "filter_by_coverage", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "filter_by_coverage", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "filter_by_coverage", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "filter_by_coverage", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "filter_by_coverage", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "filter_by_coverage", "account"),
    shell:
        """
        (
        filterbycoverage.sh \\
            in={input.fasta} \\
            cov={input.covstats} \\
            out={output.fasta} \\
            outd={output.removed_names} \\
            minc={params.minc} \\
            minp={params.minp} \\
            minr={params.minr} \\
            minl={params.minl} \\
            trim={params.trim} \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """



localrules:
    finalize_contigs,

rule finalize_contigs:
    input:
        "samples/{sample}/assembly/{sample}_final_contigs.fasta",
    output:
        "samples/{sample}/assembly/{sample}.fasta",
    log:
        "logs/samples/{sample}/assembly/finalize_contigs.log",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "account"),
    shell:
        """
        (cp {input} {output}) 1>{log} 2>&1
        """


rule calculate_contigs_stats:
    input:
        get_assembly,
    output:
        "samples/{sample}/assembly/contig_stats/final_contig_stats.txt",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "calculate_contigs_stats", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "calculate_contigs_stats", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "calculate_contigs_stats", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "calculate_contigs_stats", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "calculate_contigs_stats", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "calculate_contigs_stats", "account"),
    log:
        "logs/samples/{sample}/assembly/contig_stats/contig_stats_final.log",
    benchmark:
        "benchmarks/samples/{sample}/assembly/contig_stats/contig_stats_final.tsv",
    shell:
        """
        (stats.sh in={input} format=3 out={output} -Xmx{resources.java_mem}M) 1>{log} 2>&1
        """


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_final_contigs:
    input:
        unpack(lambda wildcards: get_quality_controlled_reads(wildcards, as_dict=True)),
        target="samples/{sample_contigs}/assembly/{sample_contigs}.fasta",
    output:
        bam=temp("samples/{sample_contigs}/sequence_alignment/{sample}.bam"),
    params:
        command = lambda wildcards, input, output, threads, resources: align_reads_command(
            wildcards, input, output, threads, resources
        ),
    benchmark:
        "benchmarks/samples/{sample_contigs}/sequence_alignment/align_reads_from_{sample}.tsv",
    log:
        "logs/samples/{sample_contigs}/sequence_alignment/align_reads_from_{sample}.log",
    conda:
        "../envs/minimap.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "mapping", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "mapping", "account"),
    shell:
        """
        ({params.command}) > {log} 2>&1
        """


rule pileup_contigs_sample:
    input:
        fasta=get_assembly,
        bam="samples/{sample}/sequence_alignment/{sample}.bam",
    output:
        covhist="samples/{sample}/assembly/contig_stats/postfilter_coverage_histogram.txt",
        covstats="samples/{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
        bincov="samples/{sample}/assembly/contig_stats/postfilter_coverage_binned.txt",
    params:
        pileup_secondary=(
            "t"
            if config["count_multi_mapped_reads"]
            else "f"
        ),
        minmapq=config["minimum_map_quality"],
    benchmark:
        "benchmarks/samples/{sample}/assembly/contig_stats/pileup_contigs_sample.tsv",
    log:
        "logs/samples/{sample}/assembly/contig_stats/pileup_contigs_sample.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "pileup", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "mem_mb"),
        java_mem        = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "java_mem"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "pileup", "account"),
    shell:
        """
        (
        pileup.sh \\
            ref={input.fasta} \\
            in={input.bam} \\
            threads={threads} \\
            -Xmx{resources.java_mem}M \\
            covstats={output.covstats} \\
            hist={output.covhist} \\
            concise=t \\
            minmapq={params.minmapq} \\
            secondary={params.pileup_secondary} \\
            bincov={output.bincov}
        ) 1>{log} 2>&1
        """


rule samtools_stats_contigs_sample:
    input:
        fasta=get_assembly,
        bam="samples/{sample}/sequence_alignment/{sample}.bam",
    output:
        stats="samples/{sample}/assembly/contig_stats/postfilter_samtools_stats.txt",
    benchmark:
        "benchmarks/samples/{sample}/assembly/contig_stats/samtools_stats_contigs_sample.tsv",
    log:
        "logs/samples/{sample}/assembly/contig_stats/samtools_stats_contigs_sample.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "samtools_stats_contigs_sample", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "samtools_stats_contigs_sample", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "samtools_stats_contigs_sample", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "samtools_stats_contigs_sample", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "samtools_stats_contigs_sample", "account"),
    shell:
        """
        (samtools stats {input.bam} 1> {output.stats}) 1>{log} 2>&1
        """


rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    log:
        "logs/{file}.index.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wildcards: get_resource(wildcards, None, 1, "create_bam_index", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "create_bam_index", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "create_bam_index", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "create_bam_index", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "create_bam_index", "account"),
    shell:
        """
        (samtools index {input}) 1>{log} 2>&1
        """


rule predict_genes:
    input:
        get_assembly,
    output:
        fna="samples/{sample}/annotation/predicted_genes/{sample}.fna",
        faa="samples/{sample}/annotation/predicted_genes/{sample}.faa",
        gff="samples/{sample}/annotation/predicted_genes/{sample}.gff",
    conda:
        "../envs/prodigal.yaml"
    log:
        "logs/samples/{sample}/annotation/predicted_genes/prodigal.log",
    benchmark:
        "benchmarks/samples/{sample}/annotation/predicted_genes/prodigal.tsv",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "predict_genes", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "predict_genes", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "predict_genes", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "predict_genes", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "predict_genes", "account"),
    shell:
        """
        (prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff) 1>{log} 2>&1
        """


localrules:
    get_contigs_from_gene_names,


rule get_contigs_from_gene_names:
    input:
        faa="samples/{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        tsv="samples/{sample}/annotation/predicted_genes/{sample}.tsv",
    log:
        "logs/samples/{sample}/annotation/predicted_genes/get_contigs_from_gene_names.log",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "account"),
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
            "samples/{sample}/assembly/contig_stats/final_contig_stats.txt", sample=SAMPLES
        ),
        gene_tables=expand(
            "samples/{sample}/annotation/predicted_genes/{sample}.tsv", sample=SAMPLES
        ),
        mapping_stats=expand(
            "samples/{sample}/assembly/contig_stats/postfilter_coverage_stats.txt", sample=SAMPLES,
        ),
        samtools_stats=expand(
            "samples/{sample}/assembly/contig_stats/postfilter_samtools_stats.txt", sample=SAMPLES
        )
    output:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    params:
        samples=SAMPLES,
    conda:
        "../envs/python.yaml"
    log:
        "logs/assembly/combine_contig_stats.log",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "account"),
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
        "logs/assembly/build_assembly_report.log",
    threads: lambda wildcards: get_resource(wildcards, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "account"),
    script:
        "../report/assembly_report.py"


