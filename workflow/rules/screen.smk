
rule generate_sketch:
    input:
        unpack(get_input_fastq),
    output:
        "Intermediate/screen/sketches/{sample}.sketch.gz",
    log:
        "logs/screen/{sample}.make_sketch.log",
    benchmark:
        "benchmarks/screen/{sample}.make_sketch.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "account"),
    shell:
        """
        (
        bbsketch.sh \\
            in={input[0]} \\
            samplerate=0.5 \\
            minkeycount=2 \\
            out={output} \\
            blacklist=nt \\
            ssu=f \\
            name0={wildcards.sample} \\
            depth=t \\
            overwrite=t \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """
        # take only one read


rule compare_sketch:
    input:
        expand(rules.generate_sketch.output, sample=SAMPLES),
    output:
        "QC/screen/sketch_comparison.tsv.gz",
    priority: 100
    log:
        "logs/screen/compare_sketch.log",
    benchmark:
        "benchmarks/screen/compare_sketch.tsv",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "mem_mb"),
        java_mem        = lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "java_mem"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "account"),
    shell:
        """
        (
        comparesketch.sh alltoall \\
            format=3 out={output} \\
            records=5000 \\
            {input} \\
            -Xmx{resources.java_mem}M
        ) 1>{log} 2>&1
        """


