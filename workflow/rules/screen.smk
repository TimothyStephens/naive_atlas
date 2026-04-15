
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
        lambda wc, input, attempt: get_all_resources(wc, input, attempt, "initialize_qc", java_mem_factor=0.85)
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
        lambda wc, input, attempt: get_all_resources(wc, input, attempt, "initialize_qc", java_mem_factor=0.85)
    shell:
        "comparesketch.sh alltoall "
        " format=3 out={output} "
        " records=5000 "
        " {input} "
        " -Xmx{resources.java_mem}M "
        " &> {log}"


#        sendsketch.sh sample2.sketch printdepth2=t level=2 printqfname=f printvolume=t color=f out
