bin_quality_input_folder = "{sample}/binning/{binner}/bins"
busco_lineages = ["archaea_odb10", "bacteria_odb10", "eukaryota_odb10"]


rule calculate_stats:
    input:
        bin_quality_input_folder,
    output:
        "{sample}/binning/{binner}/genome_stats.tsv",
    threads: config["threads"]
    params:
        extension=".fasta",
    run:
        from utils.genome_stats import get_many_genome_stats

        filenames = list(Path(input[0]).glob("*" + params.extension))

        get_many_genome_stats(filenames, output[0], threads)


rule combine_bin_stats:
    input:
        expand(
            "{sample}/binning/{{binner}}/genome_stats.tsv",
            sample=SAMPLES,
        ),
    output:
        "Binning/{binner}/raw_bins/genome_stats.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/{binner}/combine_stats.log",
    run:
        try:
            from utils.io import pandas_concat

            pandas_concat(input, output[0])

        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


#### Checkm2 ####
rule run_checkm2:
    input:
        fasta_dir=bin_quality_input_folder,
        db=rules.checkm2_download_db.output,
    output:
        table="{sample}/binning/{binner}/checkm2_report.tsv",
        faa=directory("{sample}/binning/{binner}/faa"),
    params:
        lowmem=" --lowmem " if config["mem"] < 10 else "",
        dir=lambda wc, output: Path(output.table).parent / "checkm2",
        tmpdir=lambda wc: f"{config['tmpdir']}/checkm/{wc.sample}_{wc.binner}",
    conda:
        "../envs/checkm2.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkm2.log",
    benchmark:
        "logs/benchmarks/checkm2/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        "("
        "rm -fr {params.tmpdir}; "
        "mkdir -p {params.tmpdir}; "
        "cp {CHECKM2_PSEUDO_DATA} {input.fasta_dir}/*.fasta {params.tmpdir}/; "
        "checkm2 predict"
        " --threads {threads}"
        " {params.lowmem}"
        " --force"
        " --allmodels"
        " -x .fasta"
        " --tmpdir {resources.tmpdir}"
        " --input {params.tmpdir}"
        " --output-directory {params.dir}; "
        "rm -fr {params.tmpdir}; "
        "sed -i -e '/CHECKM2_PSEUDO_DATA/d' {params.dir}/quality_report.tsv; "
        "cp {params.dir}/quality_report.tsv {output.table}; "
        "rm {params.dir}/protein_files/CHECKM2_PSEUDO_DATA.faa; "
        "mv {params.dir}/protein_files {output.faa}"
        ") 1> {log} 2>&1"


#### GUNC ####
rule run_gunc:
    input:
        db=rules.download_gunc.output[0].format(**config),
        fasta_dir="{sample}/binning/{binner}/faa",
    output:
        table="{sample}/binning/{binner}/gunc_output.tsv",
        folder=directory("{sample}/binning/{binner}/gunc"),
    params:
        extension=".faa",
    conda:
        "../envs/gunc.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/gunc.log",
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        " mkdir {output.folder} 2> {log}"
        " ;\n"
        " gunc run "
        " --threads {threads} "
        " --gene_calls "
        " --db_file {input.db} "
        " --input_dir {input.fasta_dir} "
        " --temp_dir {resources.tmpdir} "
        " --file_suffix {params.extension} "
        " --out_dir {output.folder} &>> {log} "
        " ;\n "
        " cp {output.folder}/*.tsv {output.table} 2>> {log}"


##### BUSCO  #########
rule run_busco:
    input:
        fasta_dir=bin_quality_input_folder,
        db=BUSCODIR,
    output:
        "{sample}/binning/{binner}/busco_{lineage}.tsv",
    params:
        tmpdir=lambda wc: f"{config['tmpdir']}/busco/{wc.sample}_{wc.binner}_{wc.lineage}",
    conda:
        "../envs/busco.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/busco_{lineage}.log",
    benchmark:
        "logs/benchmarks/busco/{sample}_{binner}_{lineage}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        """
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export PYTHONPATH="$CONDA_PREFIX/lib/python3.7/site-packages"
        rm -fr {params.tmpdir}
        busco -i {input.fasta_dir} \
          --lineage {wildcards.lineage} \
          -m genome \
          --out_path {params.tmpdir} \
          -o output \
          --download_path {input.db} \
          -c {threads} \
          --offline &>{log}
        awk -F'\\t' '{{ if(NR==1){{print "Bin Id\\t{wildcards.lineage}_Complete\\t{wildcards.lineage}_Single\\t{wildcards.lineage}_Duplicated\\t{wildcards.lineage}_Fragmented\\t{wildcards.lineage}_Missing"}}else{{gsub(".fasta","",$1); print $1"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7}} }}' {params.tmpdir}/output/batch_summary.txt 1>{output} 2>>{log}
        rm -fr {params.tmpdir}
        """



##### CheckV  #########
rule run_checkv:
    input:
        fasta_dir=bin_quality_input_folder,
        db=CHECKVDIR,
    output:
        "{sample}/binning/{binner}/checkv.tsv",
    params:
        tmpdir=lambda wc: f"{config['tmpdir']}/checkv/{wc.sample}_{wc.binner}",
    conda:
        "../envs/checkv.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkv.log",
    benchmark:
        "logs/benchmarks/checkv/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        """
        mkdir -p $(dirname {params.tmpdir})
        rm -fr {params.tmpdir}*
        find {input.fasta_dir} -maxdepth 1 -name "*.fasta" -exec basename {{}} \; \
          | while read FILE; do
              echo ">${{FILE%*.fasta}}"        >> {params.tmpdir}.fna
              grep -v '>' {input.fasta_dir}/${{FILE}} >> {params.tmpdir}.fna
            done &> {log}
            checkv end_to_end -t {threads} -d {input.db}/checkv-db-v1.5 {params.tmpdir}.fna {params.tmpdir}.checkv &>> {log}
            awk -F'\\t' 'BEGIN{{print "bin_id\\tcheckv_completeness\\tcheckv_contamination\\tcheckv_quality\\tcheckv_miuvig_quality"}} NR>1 {{print $1"\\t"$10"\\t"$12"\\t"$8"\\t"$9}}' \
              {params.tmpdir}.checkv/quality_summary.tsv 1>{output} 2>>{log}
        rm -fr {params.tmpdir}*
        """


##### PLASME  #########
rule run_plasme:
    input:
        fasta_dir=bin_quality_input_folder,
        db=PLASMEDONE,
    output:
        "{sample}/binning/{binner}/plasme.tsv",
    params:
        tmpdir=lambda wc: f"{config['tmpdir']}/plasme/{wc.sample}_{wc.binner}",
    conda:
        "../envs/plasme.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/plasme.log",
    benchmark:
        "logs/benchmarks/plasme/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        """
        mkdir -p $(dirname {params.tmpdir})
        rm -fr {params.tmpdir}*
        find {input.fasta_dir} -maxdepth 1 -name "*.fasta" -exec basename {{}} \; \
          | while read FILE; do
              echo ">${{FILE%*.fasta}}"        >> {params.tmpdir}.fna
              grep -v '>' {input.fasta_dir}/${{FILE}} >> {params.tmpdir}.fna
            done &> {log}
            PLASMe=$(which PLASMe.py)
            python $PLASMe --thread {threads} --database $CONDA_PREFIX/PLASMe/DB --temp {params.tmpdir}.plasme_tmp {params.tmpdir}.fna {params.tmpdir}.plasme &>> {log}
            awk -F'\\t' 'BEGIN{{print "bin_id\\tplasme_reference\\tplasme_order\\tplasme_evidence\\tplasme_score"}} NR>1 {{print $1"\\t"$3"\\t"$4"\\t"$5"\\t"$6}}' \
              {params.tmpdir}.plasme_report.csv 1>{output} 2>>{log}
        rm -fr {params.tmpdir}*
        """


## Combine bin stats
localrules:
    build_bin_report,
    combine_checkm2,
    combine_busco,
    combine_gunc,


rule combine_gunc:
    input:
        expand(
            "{sample}/binning/{{binner}}/gunc_output.tsv",
            sample=SAMPLES,
        ),
    output:
        bin_table="Binning/{binner}/raw_bins/gunc_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/{binner}/combine_gunc.log",
    run:
        try:
            from utils.io import pandas_concat

            pandas_concat(input, output[0])

        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


rule combine_checkm2:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/checkm2_report.tsv",
            sample=SAMPLES,
        ),
    output:
        bin_table="Binning/{binner}/raw_bins/checkm2_quality_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/combine_stats_{binner}_checkm2.log",
    script:
        "../scripts/combine_checkm2.py"


rule combine_busco:
    input:
        completeness_files=[ f"{x}/binning/{{binner}}/busco_{y}.tsv" for x in SAMPLES for y in busco_lineages ],
    output:
        bin_table="Binning/{binner}/raw_bins/busco_quality_report.tsv",
    params:
        samples=SAMPLES,
        lineages=busco_lineages,
    log:
        "logs/binning/combine_stats_{binner}_busco.log",
    script:
        "../scripts/combine_busco.py"


rule combine_checkv:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/checkv.tsv",
            sample=SAMPLES,
        ),
    output:
        bin_table="Binning/{binner}/raw_bins/checkv_quality_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/combine_stats_{binner}_checkv.log",
    script:
        "../scripts/combine_checkv.py"


rule combine_plasme:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/plasme.tsv",
            sample=SAMPLES,
        ),
    output:
        bin_table="Binning/{binner}/raw_bins/plasme_quality_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/combine_stats_{binner}_plasme.log",
    script:
        "../scripts/combine_plasme.py"


localrules:
    get_bin_filenames,


rule get_bin_filenames:
    input:
        dirs=expand(
            "{sample}/binning/{{binner}}/bins",
            sample=SAMPLES,
        ),
        protein_dirs=expand(
            "{sample}/binning/{{binner}}/faa",
            sample=SAMPLES,
        ),
    output:
        filenames="Binning/{binner}/raw_bins/paths.tsv",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io


        def get_list_of_files(dirs, pattern):
            fasta_files = []

            # searh for fasta files (.f*) in all bin folders
            for dir in dirs:
                dir = Path(dir)
                fasta_files += list(dir.glob(pattern))

            filenames = pd.DataFrame(fasta_files, columns=["Filename"])
            filenames.index = filenames.Filename.apply(io.simplify_path)
            filenames.index.name = "Bin"

            filenames.sort_index(inplace=True)

            return filenames


        fasta_filenames = get_list_of_files(input.dirs, "*.f*")
        faa_filenames = get_list_of_files(input.protein_dirs, "*.faa")

        assert all(
            faa_filenames.index == fasta_filenames.index
        ), "faa index and faa index are nt the same"

        faa_filenames.columns = ["Proteins"]

        filenames = pd.concat((fasta_filenames, faa_filenames), axis=1)

        filenames.to_csv(output.filenames, sep="\t")


localrules:
    all_contigs2bins,


rule all_contigs2bins:
    input:
        expand(
            "{sample}/binning/{{binner}}/cluster_attribution.tsv",
            sample=SAMPLES,
        ),
    output:
        temp("Binning/{binner}/contigs2bins.tsv.gz"),
    run:
        from utils.io import cat_files

        cat_files(input, output[0], gzip=True)


def quality_filter_bins_input(wildcards):
    "Specify input files for quality_filter_bins rule"

    input_files = dict(
        paths=rules.get_bin_filenames.output.filenames,
        stats="Binning/{binner}/raw_bins/genome_stats.tsv",
        quality="Binning/{binner}/raw_bins/checkm2_quality_report.tsv",
        busco="Binning/{binner}/raw_bins/busco_quality_report.tsv",
        checkv="Binning/{binner}/raw_bins/checkv_quality_report.tsv",
        plasme="Binning/{binner}/raw_bins/plasme_quality_report.tsv",
        gunc="Binning/{binner}/raw_bins/gunc_report.tsv",
    )

    # check if gunc is in config file
    filter_chimieric_bins = config["filter_chimieric_bins"]
    assert (
        type(filter_chimieric_bins) == bool
    ), f"filter_chimieric_bins in config file must be a boolean, got {filter_chimieric_bins}"
    if not filter_chimieric_bins:
        del input_files["gunc"]

    # replace wildcards
    for key in input_files:
        input_files[key] = input_files[key].format(binner=wildcards.binner)

    return input_files


rule quality_filter_bins:
    input:
        unpack(quality_filter_bins_input),
    output:
        info="Binning/{binner}/filtered_bin_info.tsv",
        paths="Binning/{binner}/filtered_bins_paths.txt",
    threads: 1
    log:
        "logs/binning/{binner}/filter_bins.log",
    params:
        filter_criteria=config["genome_filter_criteria"],
    script:
        "../scripts/filter_genomes.py"
