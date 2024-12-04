from glob import glob


rule get_metabat_depth_file_one_sample:
    input:
        "{sample}/sequence_alignment/{sample_reads}.bam",
    output:
        "{sample}/binning/coverage/{sample_reads}.metabat_depth.txt",
    benchmark:
        "{sample}/logs/benchmarks/binning/coverage/{sample_reads}.txt"
    log:
        "{sample}/logs/binning/coverage/{sample_reads}.log",
    conda:
        "../envs/metabat.yaml"
    threads: config["simplejob_threads"]  # multithreaded trough OMP_NUM_THREADS
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    params:
        minid=config["cobinning_readmapping_id"] * 100,
    priority: 100
    shell:
        "jgi_summarize_bam_contig_depths "
        " --percentIdentity {params.minid} "
        " --outputDepth {output} "
        " {input} &> {log} "


rule get_metabat_depth_file_combine:
    input:
        depths=lambda wc: expand(
            "{sample}/binning/coverage/{sample_reads}.metabat_depth.txt",
            sample_reads=get_alls_samples_of_group(wc),
            sample=wc.sample,
        ),
    output:
        "{sample}/binning/coverage/metabat_depth.txt",
    benchmark:
        "{sample}/logs/benchmarks/binning/coverage/metabat_depth.txt"
    log:
        "{sample}/logs/binning/coverage/metabat.log",
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    run:
        import pandas as pd
        from functools import reduce
        
        # Load each depth file and merge using columns which are identical across files: 'contigName', 'contigLen', 'totalAvgDepth'
        data_frames = [ pd.read_table(file_name, sep='\t') for file_name in input.depths ]
        df_merged = reduce(lambda left,right: pd.merge(left,right, on=['contigName', 'contigLen', 'totalAvgDepth'], how='outer'), data_frames)
        df_merged.to_csv(output, sep='\t', index=False, na_rep=0.0)


rule binning_prokaryotic:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=get_assembly,
        dbdir=rules.veba_download.output.dbdir,
    output:
        bins="{sample}/binning/veba/1_prokaryotic/{sample}/output/genomes",
        unbinned="{sample}/binning/veba/1_prokaryotic/{sample}/output/unbinned.fasta",
        stats="{sample}/binning/veba/1_prokaryotic/{sample}/output/genome_statistics.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_prokaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_prokaryotic"]["minimum_genome_length"],
        checkm2_completeness=config["veba_prokaryotic"]["checkm2_completeness"],
        checkm2_contamination=config["veba_prokaryotic"]["checkm2_contamination"],
        n_iter=config["veba_prokaryotic"]["n_iter"],
        output_path="{sample}/binning/veba/1_prokaryotic",
        output_prefix="{sample}",
    benchmark:
        "{sample}/logs/benchmarks/binning/veba/{sample}.prokaryotic.txt"
    log:
        "{sample}/logs/binning/veba/{sample}.prokaryotic.txt",
    conda:
        "../envs/VEBA-binning-prokaryotic_env.yml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["large_memory"],
        time=config["large_runtime"],
    shell:
        """
        {params.workflow_folder}/scripts/veba/binning-prokaryotic.py \
            --fasta {input.contigs} \
            --name {params.output_prefix} \
            --coverage {input.depth_file} \
            --project_directory {params.output_path} \
            --veba_database {input.dbdir} \
            --n_jobs {threads} \
            --minimum_contig_length {params.minimum_contig_length} \
            --minimum_genome_length {params.minimum_genome_length} \
            --checkm2_completeness {params.checkm2_completeness} \
            --checkm2_contamination {params.checkm2_contamination} \
            --skip_concoct \
            --n_iter {params.n_iter} \
            &> {log}
        """


rule binning_eukaryotic:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=rules.binning_prokaryotic.output.unbinned,
        dbdir=rules.veba_download.output.dbdir,
    output:
        bins="{sample}/binning/veba/2_eukaryotic/{sample}/output/genomes",
        unbinned="{sample}/binning/veba/2_eukaryotic/{sample}/output/unbinned.fasta",
        stats="{sample}/binning/veba/2_eukaryotic/{sample}/output/genome_statistics.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_eukaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_eukaryotic"]["minimum_genome_length"],
        busco_completeness=config["veba_eukaryotic"]["busco_completeness"],
        busco_contamination=config["veba_eukaryotic"]["busco_contamination"],
        output_path="{sample}/binning/veba/2_eukaryotic",
        output_prefix="{sample}",
    benchmark:
        "{sample}/logs/benchmarks/binning/veba/{sample}.eukaryotic.txt"
    log:
        "{sample}/logs/binning/veba/{sample}.eukaryotic.txt",
    conda:
        "../envs/VEBA-binning-eukaryotic_env.yml"
    threads: config["large_threads"]
    resources:
        mem=config["large_memory"],
        time=config["large_runtime"],
    shell:
        """
        {params.workflow_folder}/scripts/veba/binning-eukaryotic.py \
            --fasta {input.contigs} \
            --name {params.output_prefix} \
            --coverage {input.depth_file} \
            --project_directory {params.output_path} \
            --veba_database {input.dbdir} \
            --n_jobs {threads} \
            --minimum_contig_length {params.minimum_contig_length} \
            --minimum_genome_length {params.minimum_genome_length} \
            --busco_completeness {params.busco_completeness} \
            --busco_contamination {params.busco_contamination} \
            &> {log}
        """


rule binning_viral:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=rules.binning_eukaryotic.output.unbinned,
        dbdir=rules.veba_download.output.dbdir,
    output:
        plastid_bins="{sample}/binning/veba/3_viral/{sample}/output/filtered_plastid_bins/genomes",
        plastid_stats="{sample}/binning/veba/3_viral/{sample}/output/filtered_plastid_bins/genomad_results.filtered.tsv",
        viral_bins="{sample}/binning/veba/3_viral/{sample}/output/filtered_viral_bins/genomes",
        viral_stats="{sample}/binning/veba/3_viral/{sample}/output/filtered_viral_bins/genomad_results.filtered.tsv",
        unbinned="{sample}/binning/veba/3_viral/{sample}/output/unbinned.fasta",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_viral"]["minimum_contig_length"],
        minimum_genome_length=config["veba_viral"]["minimum_genome_length"],
        minimum_score=config["veba_viral"]["minimum_score"],
        output_path="{sample}/binning/veba/3_viral",
        output_prefix="{sample}",
    benchmark:
        "{sample}/logs/benchmarks/binning/veba/{sample}.viral.txt"
    log:
        "{sample}/logs/binning/veba/{sample}.viral.txt",
    conda:
        "../envs/VEBA-binning-viral_env.yml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["large_memory"],
        time=config["large_runtime"],
    shell:
        """
        {params.workflow_folder}/scripts/veba/binning-viral.py \
            --fasta {input.contigs} \
            --name {params.output_prefix} \
            --coverage {input.depth_file} \
            --project_directory {params.output_path} \
            --veba_database {input.dbdir} \
            --n_jobs {threads} \
            --minimum_contig_length {params.minimum_contig_length} \
            --minimum_genome_length {params.minimum_genome_length} \
            --minimum_score {params.minimum_score} \
            &> {log}
        """


## Combine binning results
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


# Combine Prokaryotic bin paths and completness stats
localrules:
    get_prokaryotic_bins,

rule get_prokaryotic_bins:
    input:
        dirs=expand(
            "{sample}/binning/veba/1_prokaryotic/{sample}/output/genomes",
            sample=SAMPLES,
        ),
        genome_stats=expand(
            "{sample}/binning/veba/1_prokaryotic/{sample}/output/genome_statistics.tsv",
            sample=SAMPLES,
        ),
        completness_stats=expand(
            "{sample}/binning/veba/1_prokaryotic/{sample}/output/checkm2_results.filtered.tsv",
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/prokaryotic.paths.tsv",
        stats="binning/raw_bins/prokaryotic.statistics.tsv",
    log:
        "logs/binning/raw_bins/get_prokaryotic_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os.path
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        faa_filenames    = get_list_of_files(input.dirs, "*.faa")
        cds_filenames    = get_list_of_files(input.dirs, "*.ffn")
        gff_filenames    = get_list_of_files(input.dirs, "*.gff")
        rRNA_filenames   = get_list_of_files(input.dirs, "*.rRNA")
        tRNA_filenames   = get_list_of_files(input.dirs, "*.tRNA")
        
        assert all(
            genome_filenames.index == faa_filenames.index
        ), "faa index does not match fa index"
        assert all(
            genome_filenames.index == cds_filenames.index
        ), "ffn index does not match fa index"
        assert all(
            genome_filenames.index == gff_filenames.index
        ), "gff index does not match fa index"
        assert all(
            genome_filenames.index == rRNA_filenames.index
        ), "rRNA index does not match fa index"
        assert all(
            genome_filenames.index == tRNA_filenames.index
        ), "tRNA index does not match fa index"
        
        faa_filenames.columns  = ["Proteins"]
        cds_filenames.columns  = ["CDS"]
        gff_filenames.columns  = ["GFF"]
        rRNA_filenames.columns = ["rRNA"]
        tRNA_filenames.columns = ["tRNA"]
        
        filenames = pd.concat((genome_filenames, faa_filenames, cds_filenames, gff_filenames, rRNA_filenames, tRNA_filenames), axis=1)
        
        filenames.to_csv(output.filenames, sep="\t")
        
        # Load each genome stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.genome_stats if os.path.isfile(file_name) ]
        df_merged_1 = pd.concat(data_frames, axis=1)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) ]
        df_merged_2 = pd.concat(data_frames, axis=1)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=False, na_rep=0.0)


## Combine Eukaryotic bin paths and completness stats
localrules:
    get_eukaryotic_bins,

rule get_eukaryotic_bins:
    input:
        dirs=expand(
            "{sample}/binning/veba/2_eukaryotic/{sample}/output/genomes",
            sample=SAMPLES,
        ),
        genome_stats=expand(
            "{sample}/binning/veba/2_eukaryotic/{sample}/output/genome_statistics.tsv",
            sample=SAMPLES,
        ),
        completness_stats=expand(
            "{sample}/binning/veba/2_eukaryotic/{sample}/output/busco_results.filtered.tsv",
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/eukaryotic.paths.tsv",
        stats="binning/raw_bins/eukaryotic.statistics.tsv",
    log:
        "logs/binning/raw_bins/get_eukaryotic_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os.path
        
        genome_filenames  = get_list_of_files(input.dirs, "*.fa")
        faa_filenames     = get_list_of_files(input.dirs, "*.faa")
        cds_filenames     = get_list_of_files(input.dirs, "*.ffn")
        gff_filenames     = get_list_of_files(input.dirs, "*.gff")
        rRNA_filenames    = get_list_of_files(input.dirs, "*.rRNA")
        tRNA_filenames    = get_list_of_files(input.dirs, "*.tRNA")
        seqType_filenames = get_list_of_files(input.dirs, "*.seq_type.tsv")
        
        assert all(
            genome_filenames.index == faa_filenames.index
        ), "faa index does not match fa index"
        assert all(
            genome_filenames.index == cds_filenames.index
        ), "ffn index does not match fa index"
        assert all(
            genome_filenames.index == gff_filenames.index
        ), "gff index does not match fa index"
        assert all(
            genome_filenames.index == rRNA_filenames.index
        ), "rRNA index does not match fa index"
        assert all(
            genome_filenames.index == tRNA_filenames.index
        ), "tRNA index does not match fa index"
        assert all(
            genome_filenames.index == seqType_filenames.index
        ), "seq_type index does not match fa index"
        
        faa_filenames.columns  = ["Proteins"]
        cds_filenames.columns  = ["CDS"]
        gff_filenames.columns  = ["GFF"]
        rRNA_filenames.columns = ["rRNA"]
        tRNA_filenames.columns = ["tRNA"]
        seqType_filenames.columns = ["seq_type"]
        
        filenames = pd.concat((genome_filenames, faa_filenames, cds_filenames, gff_filenames, rRNA_filenames, tRNA_filenames, seqType_filenames), axis=1)
        
        filenames.to_csv(output.filenames, sep="\t")
        
        # Load each genome stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.genome_stats if os.path.isfile(file_name) ]
        df_merged_1 = pd.concat(data_frames, axis=1)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0, header=[0,1]) for file_name in input.completness_stats if os.path.isfile(file_name) ]
        df_merged_2 = pd.concat(data_frames, axis=1)
        # Needed since the BUSCO results have two row header (join into one line)
        df_merged_2.columns = df_merged_2.columns.map('-'.join)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=False, na_rep=0.0)


## Combine Virus and Plasmid bin paths and completness stats
localrules:
    get_virus_bins,
    get_plasmid_bins,

rule get_virus_bins:
    input:
        dirs=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_viral_bins/genomes",
            sample=SAMPLES,
        ),
        genome_stats=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_viral_bins/genome_statistics.tsv",
            sample=SAMPLES,
        ),
        completness_stats=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_viral_bins/genomad_results.filtered.tsv",
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/viral.paths.tsv",
        stats="binning/raw_bins/viral.genome_statistics.tsv",
    log:
        "logs/binning/raw_bins/get_virus_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os.path
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        faa_filenames    = get_list_of_files(input.dirs, "*.faa")
        cds_filenames    = get_list_of_files(input.dirs, "*.ffn")
        gff_filenames    = get_list_of_files(input.dirs, "*.gff")
        
        assert all(
            genome_filenames.index == faa_filenames.index
        ), "faa index does not match fa index"
        assert all(
            genome_filenames.index == cds_filenames.index
        ), "ffn index does not match fa index"
        assert all(
            genome_filenames.index == gff_filenames.index
        ), "gff index does not match fa index"
        
        faa_filenames.columns  = ["Proteins"]
        cds_filenames.columns  = ["CDS"]
        gff_filenames.columns  = ["GFF"]
        
        filenames = pd.concat((genome_filenames, faa_filenames, cds_filenames, gff_filenames), axis=1)
        
        filenames.to_csv(output.filenames, sep="\t")
        
        # Load each genome stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.genome_stats if os.path.isfile(file_name) ]
        df_merged_1 = pd.concat(data_frames, axis=1)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) ]
        df_merged_2 = pd.concat(data_frames, axis=1)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=False, na_rep=0.0)


rule get_plasmid_bins:
    input:
        dirs=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_plasmid_bins/genomes",
            sample=SAMPLES,
        ),
        genome_stats=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_plasmid_bins/genome_statistics.tsv",
            sample=SAMPLES,
        ),
        completness_stats=expand(
            "{sample}/binning/veba/3_viral/{sample}/output/filtered_plasmid_bins/genomad_results.filtered.tsv",
            sample=SAMPLES,
        ),
    output:
        filenames="binning/raw_bins/plasmid.paths.tsv",
        stats="binning/raw_bins/plasmid.genome_statistics.tsv",
    log:
        "logs/binning/raw_bins/get_plasmid_bins.log",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io
        import os.path
        
        genome_filenames = get_list_of_files(input.dirs, "*.fa")
        faa_filenames    = get_list_of_files(input.dirs, "*.faa")
        cds_filenames    = get_list_of_files(input.dirs, "*.ffn")
        gff_filenames    = get_list_of_files(input.dirs, "*.gff")
        
        assert all(
            genome_filenames.index == faa_filenames.index
        ), "faa index does not match fa index"
        assert all(
            genome_filenames.index == cds_filenames.index
        ), "ffn index does not match fa index"
        assert all(
            genome_filenames.index == gff_filenames.index
        ), "gff index does not match fa index"
        
        faa_filenames.columns  = ["Proteins"]
        cds_filenames.columns  = ["CDS"]
        gff_filenames.columns  = ["GFF"]
        
        filenames = pd.concat((genome_filenames, faa_filenames, cds_filenames, gff_filenames), axis=1)
        
        filenames.to_csv(output.filenames, sep="\t")
        
        # Load each genome stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.genome_stats if os.path.isfile(file_name) ]
        df_merged_1 = pd.concat(data_frames, axis=1)
        
        # Load each completness stats file and concat
        data_frames = [ pd.read_table(file_name, sep='\t', index_col=0) for file_name in input.completness_stats if os.path.isfile(file_name) ]
        df_merged_2 = pd.concat(data_frames, axis=1)
        
        df_merged = pd.merge(df_merged_1, df_merged_2, how='outer', left_index=True, right_index=True)
        df_merged.to_csv(output.stats, sep='\t', index=False, na_rep=0.0)


