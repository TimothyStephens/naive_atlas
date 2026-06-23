


#################################
####                         ####
#### Combine binning results ####
####                         ####
#################################

def get_list_of_files(dirs, pattern):
    from utils import io
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


# Combine Prokaryotic bin paths and completness stats
localrules:
    get_prokaryotic_bins,

rule get_prokaryotic_bins:
    input:
        dirs=expand(
            rules.binning_prokaryotic_checkm2.output.bins,
            sample=SAMPLES,
        ),
        genome_stats=expand(
            rules.binning_prokaryotic_genome_stats.output.stats,
            sample=SAMPLES,
        ),
        completness_stats=expand(
            rules.binning_prokaryotic_checkm2.output.quality,
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/prokaryotic.paths.tsv",
        genome_names="binning/raw_bins/prokaryotic.genome.paths.tsv",
        stats="binning/raw_bins/prokaryotic.statistics.tsv",
    params:
        dir="binning/raw_bins",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/get_prokaryotic_bins.log",
    conda:
        "../envs/python.yaml"
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
        os.makedirs(params.dir, exist_ok=True)
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        if genome_filenames.empty:
            print("[WARNING] No Bacterial bins found!")
            Path(output.filenames).touch()
            Path(output.genome_names).touch()
            Path(output.stats).touch()
            return
        
        genome_filenames.columns  = ["Genome"]
        filenames = genome_filenames
        filenames.to_csv(output.filenames, sep="\t")
        filenames['Genome'].to_csv(output.genome_names, index=False, header=False)
        
        # Load each genome stats file and concat
        data_frames = []
        for file_name in input.genome_stats:
            if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
                t = pd.read_table(file_name, sep='\t', index_col=0)
                t['Sample'] = file_name.split(os.sep)[0]
                data_frames.append(t)
        df_merged_1 = pd.concat(data_frames, axis=0)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) and os.stat(file_name).st_size != 0 ]
        df_merged_2 = pd.concat(data_frames, axis=0)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=True, na_rep=0.0)


## Combine Eukaryotic bin paths and completness stats
localrules:
    get_eukaryotic_bins,

rule get_eukaryotic_bins:
    input:
        dirs=expand(
            rules.binning_eukaryotic_filter.output.bins,
            sample=SAMPLES,
        ),
        genome_stats=expand(
            rules.binning_eukaryotic_genome_stats.output.stats,
            sample=SAMPLES,
        ),
        completness_stats=expand(
            rules.binning_eukaryotic_filter.output.quality,
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/eukaryotic.paths.tsv",
        genome_names="binning/raw_bins/eukaryotic.genome.paths.tsv",
        stats="binning/raw_bins/eukaryotic.statistics.tsv",
    params:
        dir="binning/raw_bins",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/get_eukaryotic_bins.log",
    conda:
        "../envs/python.yaml"
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
        os.makedirs(params.dir, exist_ok=True)
        
        genome_filenames  = get_list_of_files(input.dirs, "*.fa")
        if genome_filenames.empty:
            print("[WARNING] No Eukaryotic bins found!")
            Path(output.filenames).touch()
            Path(output.genome_names).touch()
            Path(output.stats).touch()
            return
        
        genome_filenames.columns  = ["Genome"]
        filenames = genome_filenames
        filenames.to_csv(output.filenames, sep="\t")
        filenames['Genome'].to_csv(output.genome_names, index=False, header=False)
        
        # Load each genome stats file and concat
        data_frames = []
        for file_name in input.genome_stats:
            if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
                t = pd.read_table(file_name, sep='\t', index_col=0)
                t['Sample'] = file_name.split(os.sep)[0]
                data_frames.append(t)
        df_merged_1 = pd.concat(data_frames, axis=0)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0, header=[0,1]) for file_name in input.completness_stats if os.path.isfile(file_name) and os.stat(file_name).st_size != 0 ]
        df_merged_2 = pd.concat(data_frames, axis=0)
        # Needed since the BUSCO results have two row header (join into one line)
        df_merged_2.columns = df_merged_2.columns.map('-'.join)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=True, na_rep=0.0)


## Combine Virus and Plasmid bin paths and completness stats
localrules:
    get_viral_bins,
    get_plasmid_bins,

rule get_viral_bins:
    input:
        dirs=expand(
            rules.binning_viral_filter.output.viral_bins,
            sample=SAMPLES,
        ),
        genome_stats=expand(
            rules.binning_viral_genome_stats.output.viral_stats,
            sample=SAMPLES,
        ),
        completness_stats=expand(
            rules.binning_viral_filter.output.viral_completness,
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/viral.paths.tsv",
        genome_names="binning/raw_bins/viral.genome.paths.tsv",
        stats="binning/raw_bins/viral.statistics.tsv",
    params:
        dir="binning/raw_bins",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/get_viral_bins.log",
    conda:
        "../envs/python.yaml"
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
        os.makedirs(params.dir, exist_ok=True)
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        if genome_filenames.empty:
            print("[WARNING] No Viral bins found!")
            Path(output.filenames).touch()
            Path(output.genome_names).touch()
            Path(output.stats).touch()
            return
        
        genome_filenames.columns  = ["Genome"]
        filenames = genome_filenames
        filenames.to_csv(output.filenames, sep="\t")
        filenames['Genome'].to_csv(output.genome_names, index=False, header=False)
        
        # Load each genome stats file and concat
        data_frames = []
        for file_name in input.genome_stats:
            if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
                t = pd.read_table(file_name, sep='\t', index_col=0)
                t['Sample'] = file_name.split(os.sep)[0]
                data_frames.append(t)
        df_merged_1 = pd.concat(data_frames, axis=0)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) and os.stat(file_name).st_size != 0 ]
        df_merged_2 = pd.concat(data_frames, axis=0)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=True, na_rep=0.0)


rule get_plasmid_bins:
    input:
        dirs=expand(
            rules.binning_viral_filter.output.plasmid_bins,
            sample=SAMPLES,
        ),
        genome_stats=expand(
            rules.binning_viral_genome_stats.output.plasmid_stats,
            sample=SAMPLES,
        ),
        completness_stats=expand(
            rules.binning_viral_filter.output.plasmid_completness,
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/plasmid.paths.tsv",
        genome_names="binning/raw_bins/plasmid.genome.paths.tsv",
        stats="binning/raw_bins/plasmid.statistics.tsv",
    params:
        dir="binning/raw_bins",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/get_plasmid_bins.log",
    conda:
        "../envs/python.yaml"
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
        os.makedirs(params.dir, exist_ok=True)
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        if genome_filenames.empty:
            print("[WARNING] No Plasmid bins found!")
            Path(output.filenames).touch()
            Path(output.genome_names).touch()
            Path(output.stats).touch()
            return
        
        genome_filenames.columns  = ["Genome"]
        filenames = genome_filenames
        filenames.to_csv(output.filenames, sep="\t")
        filenames['Genome'].to_csv(output.genome_names, index=False, header=False)
        
        # Load each genome stats file and concat
        data_frames = []
        for file_name in input.genome_stats:
            if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
                t = pd.read_table(file_name, sep='\t', index_col=0)
                t['Sample'] = file_name.split(os.sep)[0]
                data_frames.append(t)
        df_merged_1 = pd.concat(data_frames, axis=0)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) and os.stat(file_name).st_size != 0 ]
        df_merged_2 = pd.concat(data_frames, axis=0)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=True, na_rep=0.0)


localrules:
    get_all,

checkpoint get_all:
    input:
        paths=expand("binning/raw_bins/{lineage}.genome.paths.tsv", lineage=['prokaryotic', 'eukaryotic', 'viral', 'plasmid']),
    log:
        "logs/binning/raw_bins/get_all.log",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    output:
        touch("binning/raw_bins/all.done"),





#################################
####                         ####
####       DeReplicate       ####
####                         ####
#################################

rule run_skani:
    input:
        all_done="binning/raw_bins/all.done",
        paths="binning/raw_bins/{lineage}.genome.paths.tsv",
    output:
        "binning/raw_bins/{lineage}.distance_matrix.txt",
    log:
        "logs/binning/raw_bins/{lineage}.skani_calculation.log",
    benchmark:
        "benchmarks/binning/raw_bins/{lineage}.skani_calculation.tsv",
    threads: lambda wc: get_resource(wc, None, 1, "run_skani", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_skani", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_skani", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_skani", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_skani", "account"),
    params:
        lineage="{lineage}",
        min_af=config["genome_dereplication"]["overlap"] * 100,
        extra="",
    conda:
        "../envs/skani.yaml"
    shell:
        """
        (
        sensitivity="--medium"
        if [ "{params.lineage}" == "viral" ] || [ "{params.lineage}" == "plasmid" ]; then
          sensitivity="--slow"
        fi
        
        skani triangle \\
          {params.extra} \\
          -l {input.paths} \\
          -o {output} \\
          -t {threads} \\
          --sparse --ci \\
          --min-af {params.min_af} \\
          $sensitivity
        ) 1>{log} 2>&1
        """


rule skani_2_parquet:
    input:
        rules.run_skani.output,
    output:
        "binning/raw_bins/{lineage}.genome_similarities.parquet",
    threads: lambda wc: get_resource(wc, None, 1, "skani_2_parquet", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "skani_2_parquet", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "skani_2_parquet", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "skani_2_parquet", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "skani_2_parquet", "account"),
    log:
        "logs/binning/raw_bins/{lineage}.skani_2_parquet.log",
    benchmark:
        "benchmarks/binning/raw_bins/{lineage}.skani_2_parquet.tsv",
    conda:
        "../envs/python.yaml"
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
        bin_info="binning/raw_bins/{lineage}.statistics.tsv",
    params:
        linkage_method="average",
        pre_cluster_threshold=0.925,
        threshold=config["genome_dereplication"]["ANI"],
        script="../scripts/cluster_{lineage}_species.py"
    conda:
        "../envs/species_clustering.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "cluster_species", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "cluster_species", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "cluster_species", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "cluster_species", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "cluster_species", "account"),
    log:
        "logs/binning/raw_bins/{lineage}.species_clustering.log",
    benchmark:
        "benchmarks/binning/raw_bins/{lineage}.species_clustering.tsv",
    output:
        bin_info="binning/{lineage}.bin_info.tsv",
        bins2species="binning/{lineage}.bins2species.tsv",
    script:
        "{params.script}"


rule build_bin_report:
    input:
        bin_info="binning/{lineage}.bin_info.tsv",
        bins2species="binning/{lineage}.bins2species.tsv",
    output:
        report="reports/bin_report_{lineage}.html",
    params:
        script="../report/bin_report_{lineage}.py"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_{lineage}.log",
    script:
        "{params.script}"


rule run_cdhit:
    input:
        expand("samples/{sample}/binning/veba/3_viral/2_genomad/unbinned.fasta",
            sample=SAMPLES
        ),
    output:
        "binning/raw_unbinned/combined.cdhit_est",
    log:
        "logs/binning/raw_unbinned/run_cdhit.log",
    benchmark:
        "benchmarks/binning/raw_unbinned/run_cdhit.tsv",
    threads: lambda wc: get_resource(wc, None, 1, "run_cdhit", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_cdhit", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_cdhit", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_cdhit", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "run_cdhit", "account"),
    params:
        combined="binning/raw_unbinned/combined.fa",
        prefix="binning/raw_unbinned/combined.cdhit_est",
        extra=config["unbinned_dereplication"]["cdhitest_params"],
    container:
        "docker://chrishah/cdhit:v4.8.1"
    shell:
        """
        (
        cat {input} > {params.combined}
        cd-hit-est -i {params.combined} -o {params.prefix} {params.extra} -T {threads}
        ) 1>{log} 2>&1
        """





#################################
####                         ####
#### Rename and move genomes ####
####                         ####
#################################

localrules:
    rename_genomes,
    rename_unbinned,
    move_genomes,
    move_unbinned,


rule rename_genomes:
    input:
        paths="binning/raw_bins/{lineage}.paths.tsv",
        mapping_file="binning/{lineage}.bins2species.tsv",
        genome_info="binning/{lineage}.bin_info.tsv",
    output:
        dir=directory("tmp/genomes/{lineage}"),
        mapfile_c2g="genomes/clustering/{lineage}.contig2genome.tsv",
        mapfile_g2c="genomes/clustering/{lineage}.genome2contig.tsv",
        mapfile_old2mag="genomes/clustering/{lineage}.old2newID.tsv",
        mapfile_allbins2mag="genomes/clustering/{lineage}.allbins2genome.tsv",
        genome_info="tmp/genomes/MAG_{lineage}.genome_quality.tsv",
    params:
        rename_contigs=config["rename_mags_contigs"],
        prefix="MAG_{lineage}_",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/{lineage}.rename_genomes.log",
    script:
        "../scripts/rename_genomes.py"


rule rename_unbinned:
    input:
        fa=rules.run_cdhit.output,
    output:
        fa="tmp/unbinned/Unbinned.fa",
        mapfile_c2g="genomes/clustering/unbinned.contig2genome.tsv",
        mapfile_g2c="genomes/clustering/unbinned.genome2contig.tsv",
    params:
        rename_contigs=config["rename_mags_contigs"],
        prefix="Unbinned",
        outdir="tmp/unbinned",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/rename_unbinned.log",
    script:
        "../scripts/rename_unbinned.py"


def get_binned_lineages(wildcards):
    checkpoints.get_all.get(**wildcards)
    binned_lineages = []
    for lineage in ['prokaryotic', 'eukaryotic', 'viral', 'plasmid']:
        file_name = f"binning/raw_bins/{lineage}.genome.paths.tsv"
        if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
            binned_lineages.append(lineage)
    return binned_lineages


def get_c2g(wildcards):
    return(expand("genomes/clustering/{lineage}.contig2genome.tsv",
                lineage=get_binned_lineages(wildcards)
            )
    )
def get_g2c(wildcards):
    return(expand("genomes/clustering/{lineage}.genome2contig.tsv",
                lineage=get_binned_lineages(wildcards)
            )
    )

rule combine_name_mappings:
    input:
        c2g_bins=get_c2g,
	g2c_bins=get_g2c,
        c2g_unbinned="genomes/clustering/unbinned.contig2genome.tsv",
        g2c_unbinned="genomes/clustering/unbinned.genome2contig.tsv",
    output:
        c2g_all="genomes/clustering/all.contig2genome.tsv",
        g2c_all="genomes/clustering/all.genome2contig.tsv",
        c2g_MAGs="genomes/clustering/mags.contig2genome.tsv",
        g2c_MAGs="genomes/clustering/mags.genome2contig.tsv",
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    log:
        "logs/genomes/clustering/combine_name_mappings.log",
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    shell:
        """
        cat {input.c2g_bins} {input.c2g_unbinned} > {output.c2g_all}
        cat {input.g2c_bins} {input.g2c_unbinned} > {output.g2c_all}
	# "cat"ing /dev/null prevents the rule from hanging if we dont have any MAGs. Will not affect output.
        cat {input.c2g_bins} /dev/null > {output.c2g_MAGs}
        cat {input.g2c_bins} /dev/null > {output.g2c_MAGs}
        """


def get_genome_to_move(wildcards):
    return(expand("tmp/genomes/{lineage}", 
                lineage=get_binned_lineages(wildcards)
            )
    )

checkpoint move_genomes:
    input:
        all_done="binning/raw_bins/all.done",
        dirs=get_genome_to_move,
    output:
        dir=directory("genomes/genomes"),
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/move_mags.log",
    script:
        "../scripts/move_genomes.sh"


checkpoint move_unbinned:
    input:
        fa=rules.rename_unbinned.output.fa,
    output:
        dir=directory("genomes/unbinned"),
        fa="genomes/unbinned/Unbinned.fa",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb          = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        runtime         = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "time_min"),
        slurm_partition = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "partition"),
        slurm_account   = lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "account"),
    log:
        "logs/binning/raw_bins/move_unbinned.log",
    script:
        "../scripts/move_unbinned.sh"


