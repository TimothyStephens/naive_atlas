
rule generate_sketch:
    input:
        unpack(get_input_fastq),
    output:
        "Intermediate/screen/sketches/{sample}.sketch.gz",
    log:
        "logs/screen/make_sketch/{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "generate_sketch", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "generate_sketch", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "generate_sketch", "account"),
        java_mem=lambda wc, input, attempt: int(get_resource(wc, input, attempt, "generate_sketch", "mem_mb") * 0.85),
    shell:
        "bbsketch.sh "
        "in={input[0]}"
        " samplerate=0.5"
        " minkeycount=2 "
        " out={output} "
        " blacklist=nt ssu=f name0={wildcards.sample} depth=t overwrite=t "
        " -Xmx{resources.java_mem}M "
        " &> {log}"
        # take only one read


rule compare_sketch:
    input:
        expand(rules.generate_sketch.output, sample=SAMPLES),
    output:
        "QC/screen/sketch_comparison.tsv.gz",
    priority: 100
    log:
        "logs/screen/compare_sketch.log",
    conda:
        "../envs/required_packages.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "initialize_qc", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "compare_sketch", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "compare_sketch", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "compare_sketch", "account"),
        java_mem=lambda wc, input, attempt: int(get_resource(wc, input, attempt, "compare_sketch", "mem_mb") * 0.85),
    shell:
        "comparesketch.sh alltoall "
        " format=3 out={output} "
        " records=5000 "
        " {input} "
        " -Xmx{resources.java_mem}M "
        " &> {log}"


#        sendsketch.sh sample2.sketch printdepth2=t level=2 printqfname=f printvolume=t color=f out
