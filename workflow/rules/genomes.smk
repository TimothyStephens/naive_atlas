


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
        filenames="Binning/raw_bins/prokaryotic.paths.tsv",
        genome_names="Binning/raw_bins/prokaryotic.genome.paths.tsv",
        stats="Binning/raw_bins/prokaryotic.statistics.tsv",
    log:
        "logs/Binning/raw_bins/get_prokaryotic_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
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
        filenames="Binning/raw_bins/eukaryotic.paths.tsv",
        genome_names="Binning/raw_bins/eukaryotic.genome.paths.tsv",
        stats="Binning/raw_bins/eukaryotic.statistics.tsv",
    log:
        "logs/Binning/raw_bins/get_eukaryotic_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
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
        filenames="Binning/raw_bins/viral.paths.tsv",
        genome_names="Binning/raw_bins/viral.genome.paths.tsv",
        stats="Binning/raw_bins/viral.statistics.tsv",
    log:
        "logs/Binning/raw_bins/get_viral_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
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
        filenames="Binning/raw_bins/plasmid.paths.tsv",
        genome_names="Binning/raw_bins/plasmid.genome.paths.tsv",
        stats="Binning/raw_bins/plasmid.statistics.tsv",
    log:
        "logs/Binning/raw_bins/get_plasmid_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os, os.path
        
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


checkpoint get_all:
    input:
        paths=expand("Binning/raw_bins/{lineage}.genome.paths.tsv", lineage=['prokaryotic', 'eukaryotic', 'viral', 'plasmid']),
    output:
        touch("Binning/raw_bins/all.done"),





#################################
####                         ####
####       DeReplicate       ####
####                         ####
#################################

rule run_skani:
    input:
        all_done="Binning/raw_bins/all.done",
        paths="Binning/raw_bins/{lineage}.genome.paths.tsv",
    output:
        "Binning/raw_bins/{lineage}.distance_matrix.txt",
    log:
        "logs/Binning/dereplication/{lineage}.skani_calculation.log",
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
        "Binning/raw_bins/{lineage}.genome_similarities.parquet",
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    log:
        "logs/Binning/dereplication/{lineage}.skani_2_parquet.log",
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
        dist="Binning/raw_bins/{lineage}.genome_similarities.parquet",
        bin_info="Binning/raw_bins/{lineage}.statistics.tsv",
    params:
        linkage_method="average",
        pre_cluster_threshold=0.925,
        threshold=config["genome_dereplication"]["ANI"],
        script="../scripts/cluster_{lineage}_species.py"
    conda:
        "../envs/species_clustering.yaml"
    log:
        "logs/Binning/dereplication/{lineage}.species_clustering.log",
    output:
        bin_info="Binning/{lineage}.bin_info.tsv",
        bins2species="Binning/{lineage}.bins2species.tsv",
    script:
        "{params.script}"


rule build_bin_report:
    input:
        bin_info="Binning/{lineage}.bin_info.tsv",
        bins2species="Binning/{lineage}.bins2species.tsv",
    output:
        report="reports/bin_report_{lineage}.html",
    params:
        script="../report/bin_report_{lineage}.py"
    conda:
        "../envs/report.yaml"
    log:
        "logs/Binning/report_{lineage}.log",
    script:
        "{params.script}"





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
        paths="Binning/raw_bins/{lineage}.paths.tsv",
        mapping_file="Binning/{lineage}.bins2species.tsv",
        genome_info="Binning/{lineage}.bin_info.tsv",
    output:
        dir=directory("tmp/genomes/{lineage}"),
        mapfile_contigs="genomes/clustering/{lineage}.contig2genome.tsv",
        mapfile_old2mag="genomes/clustering/{lineage}.old2newID.tsv",
        mapfile_allbins2mag="genomes/clustering/{lineage}.allbins2genome.tsv",
        genome_info="tmp/genomes/MAG_{lineage}.genome_quality.tsv",
    params:
        rename_contigs=config["rename_mags_contigs"],
        prefix="MAG_{lineage}_",
    shadow:
        "shallow"
    log:
        "logs/genomes/clustering/{lineage}.rename_genomes.log",
    script:
        "../scripts/rename_genomes.py"


rule rename_unbinned:
    input:
        unbinned="{sample}/binning/veba/3_viral/2_genomad/unbinned.fasta",
    output:
        dir=directory("tmp/unbinned/{sample}"),
    params:
        rename_contigs=config["rename_mags_contigs"],
        prefix="Unbinned_{sample}",
    shadow:
        "shallow"
    log:
        "logs/genomes/clustering/{sample}.rename_unbinned.log",
    script:
        "../scripts/rename_unbinned.py"


def get_binned_lineages():
    binned_lineages = []
    for lineage in ['prokaryotic', 'eukaryotic', 'viral', 'plasmid']:
        file_name = f"Binning/raw_bins/{lineage}.genome.paths.tsv"
        if os.path.isfile(file_name) and os.stat(file_name).st_size != 0:
            binned_lineages.append(lineage)
    return binned_lineages

def get_genome_to_move(wildcards):
    return(expand("tmp/genomes/{lineage}", 
                lineage=get_binned_lineages()
            )
    )

rule move_genomes:
    input:
        all_done="Binning/raw_bins/all.done",
        dirs=get_genome_to_move,
    output:
        dir=directory(GENOME_DIR),
    log:
        "logs/genomes/move_mags.log",
    script:
        "../scripts/move_genomes.sh"


rule move_unbinned:
    input:
        dirs=expand("tmp/unbinned/{sample}",
            sample=SAMPLES
        ),
    output:
        dir=directory(UNBINNED_DIR),
    log:
        "logs/genomes/move_unbinned.log",
    script:
        "../scripts/move_unbinned.sh"



