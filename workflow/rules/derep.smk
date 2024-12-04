binned_lineages = ['prokaryotic', 'eukaryotic', 'virus', 'plasmid']

rule run_skani:
    input:
        paths="binning/raw_bins/{lineage}.paths.tsv",
    output:
        "binning/raw_bins/{lineage}.distance_matrix.txt",
    log:
        "logs/binning/dereplication/{lineage}.skani_calculation.log",
    resources:
        mem_mb=config["simplejob_memory"] * 1000,
        time_min=60 * config["simplejob_runtime"],
    params:
        #preset= "medium", # fast, medium or slow
        min_af=config["genome_dereplication"]["overlap"] * 100,
        extra="",
    threads: config["simplejob_threads"]
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle "
        " {params.extra} "
        " -l {input.paths} "
        " -o {output} "
        " -t {threads} "
        " --sparse --ci "
        " --min-af {params.min_af} "
        " &> {log} "


rule skani_2_parquet:
    input:
        rules.run_skani.output,
    output:
        "binning/raw_bins/{lineage}.genome_similarities.parquet",
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    log:
        "logs/binning/dereplication/{lineage}.skani_2_parquet.log",
    run:
        try:
            skani_column_dtypes = {
                "Ref_file": "category",
                "Query_file": "category",
                "ANI": float,
                "Align_fraction_ref": float,
                "Align_fraction_query": float,
                "ANI_5_percentile": float,
                "ANI_95_percentile": float,
            }  # Ref_name        Query_name

            import pandas as pd
            from utils.io import simplify_path
            
            df = pd.read_table(input[0])
            df = pd.read_table(
                input[0],
                usecols=list(skani_column_dtypes.keys()),
                dtype=skani_column_dtypes,
            )
            
            df["Ref"] = df.Ref_file.cat.rename_categories(simplify_path)
            df["Query"] = df.Query_file.cat.rename_categories(simplify_path)
            df.to_parquet(output[0])
        
        except Exception as e:
            import traceback
            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)
            raise e


rule cluster_species:
    input:
        dist="binning/raw_bins/{lineage}.genome_similarities.parquet",
        bin_info="binning/raw_bins/{lineage}.genome_statistics.tsv",
    params:
        linkage_method="average",
        pre_cluster_threshold=0.925,
        threshold=config["genome_dereplication"]["ANI"],
    conda:
        "../envs/species_clustering.yaml"
    log:
        "logs/binning/dereplication/{lineage}.species_clustering.log",
    output:
        bin_info="binning/{lineage}.bin_info.tsv",
        bins2species="binning/{lineage}.bins2species.tsv",
    script:
        "../scripts/cluster_{lineage}_species.py"


rule build_bin_report:
    input:
        bin_info="binning/{lineage}.bin_info.tsv",
        bins2species="binning/{lineage}.bins2species.tsv",
    output:
        report="reports/bin_report_{lineage}.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_{lineage}.log",
    script:
        "../report/bin_report_{lineage}.py"

