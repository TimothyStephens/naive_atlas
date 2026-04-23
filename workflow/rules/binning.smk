from glob import glob



###############################
####                       ####
#### Prep Data for Binning ####
####                       ####
###############################
rule get_metabat_depth_file_one_sample:
    input:
        "samples/{sample}/sequence_alignment/{sample_reads}.bam",
    output:
        "samples/{sample}/binning/coverage/{sample_reads}.metabat_depth.txt",
    benchmark:
        "benchmarks/samples/{sample}/binning/coverage/{sample_reads}.txt"
    log:
        "logs/samples/{sample}/binning/coverage/{sample_reads}.log",
    conda:
        "../envs/metabat2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
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
            "samples/{{sample}}/binning/coverage/{sample_reads}.metabat_depth.txt",
            sample_reads=get_alls_samples_of_group(wc),
        ),
    output:
        depth="samples/{sample}/binning/coverage/metabat_depth.txt",
    params:
        workflow_folder=f"{workflow_folder}",
    benchmark:
        "benchmarks/samples/{sample}/binning/coverage/metabat_depth.txt"
    log:
        "logs/samples/{sample}/binning/coverage/metabat_depth.log",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        {params.workflow_folder}/scripts/veba/metabat2_coverage_combined_file.py \
            -i {input.depths} \
            -o {output.depth} \
          > {log} 2>&1
        """


rule get_maxbin_depth_file:
    input:
        depth="samples/{sample}/binning/coverage/metabat_depth.txt",
    output:
        depth="samples/{sample}/binning/coverage/maxbin_depth.txt",
    params:
        workflow_folder=f"{workflow_folder}",
    benchmark:
        "benchmarks/samples/{sample}/binning/coverage/maxbin_depth.txt"
    log:
        "logs/samples/{sample}/binning/coverage/maxbin_depth.log",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        {params.workflow_folder}/scripts/veba/maxbin_abundance_from_metabat2_coverage_file.py \
            -i {input.depth} \
            -o {output.depth} \
          > {log} 2>&1
        """





###############################
####                       ####
####  Binning Prokaryotes  ####
####                       ####
###############################

rule binning_prokaryotic_metabat:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=get_assembly,
    output:
        s2b="samples/{sample}/binning/veba/1_prokaryotic/1_metabat2/scaffolds_to_bins.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_prokaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_prokaryotic"]["minimum_genome_length"],
        output_path="samples/{sample}/binning/veba/1_prokaryotic/1_metabat2",
        output_prefix="{sample}__METABAT2",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/1_metabat2.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/1_metabat2.log",
    conda:
        "../envs/metabat2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        metabat2 \
            -i {input.contigs} \
            -o {params.output_path}/bins/bin \
            -a {input.depth_file} \
            -m {params.minimum_contig_length} \
            --minClsSize {params.minimum_genome_length} \
            -t {threads} \
            --seed 1 \
            --verbose \
        
        {params.workflow_folder}/scripts/veba/scaffolds_to_bins.py \
            -x fa \
            -i {params.output_path}/bins \
            --bin_prefix {params.output_prefix} \
            > {output.s2b}
        ) &> {log}
        """


rule binning_prokaryotic_maxbin_107:
    input:
        depth_file=rules.get_maxbin_depth_file.output,
        contigs=get_assembly,
    output:
        s2b="samples/{sample}/binning/veba/1_prokaryotic/2_maxbin2_107/scaffolds_to_bins.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_prokaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_prokaryotic"]["minimum_genome_length"],
        output_path="samples/{sample}/binning/veba/1_prokaryotic/2_maxbin2_107",
        output_prefix="{sample}__MAXBIN2-107",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/2_maxbin2_107.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/2_maxbin2_107.log",
    #conda:
    #    "../envs/maxbin2.yaml"
    container:
        "docker://timothystephens/maxbin2:2.2.7-TGSv5",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    retries: 5
    shell:
        """
        (
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        run_MaxBin.pl \
            -contig {input.contigs} \
            -out {params.output_path}/bin \
            -abund_list {input.depth_file} \
            -min_contig_length {params.minimum_contig_length} \
            -markerset 107 \
            -thread {threads} #-verbose
        
        mkdir -p {params.output_path}/bins
        
        if grep -q 'Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.' "{params.output_path}/bin.log" \
        || grep -q 'This suggests that the dataset cannot be binned (likely too few and/or small contigs), rather then it actually being an error.' "{params.output_path}/bin.log" \
        || grep -q 'Yielded 0 bins for contig (scaffold) file' "{params.output_path}/bin.log";
        then
            echo "[WARNING]  - Looks like MaxBin2-107 didnt found any prokaryotic bins, this is not a problem and expected for some samples."
            touch "{output.s2b}"
            exit 0
        fi
        
        for FP in {params.output_path}/bin.*.fasta;
        do
            GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\\n" | wc -m)
            echo "[GENOME SIZE] ${{FP}} = ${{GENOME_SIZE}}"
            
            if (( {params.minimum_genome_length} > ${{GENOME_SIZE}} ));
            then
                echo "[COPYING] ${{FP}}"
                ID_GENOME=$(basename ${{FP}} .fasta)
                mv $FP {params.output_path}/bins/${{ID_GENOME}}.fa
            fi
        done
        
        {params.workflow_folder}/scripts/veba/scaffolds_to_bins.py \
            -x fa \
            -i {params.output_path}/bins \
            --bin_prefix {params.output_prefix} \
            > {output.s2b}
        ) &> {log}
        """


rule binning_prokaryotic_maxbin_40:
    input:
        depth_file=rules.get_maxbin_depth_file.output,
        contigs=get_assembly,
    output:
        s2b="samples/{sample}/binning/veba/1_prokaryotic/3_maxbin2_40/scaffolds_to_bins.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_prokaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_prokaryotic"]["minimum_genome_length"],
        output_path="samples/{sample}/binning/veba/1_prokaryotic/3_maxbin2_40",
        output_prefix="{sample}__MAXBIN2-40",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/3_maxbin2_40.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/3_maxbin2_40.log",
    #conda:
    #    "../envs/maxbin2.yaml"
    container:
        "docker://timothystephens/maxbin2:2.2.7-TGSv5",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    retries: 5
    shell:
        """
        (
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        run_MaxBin.pl \
            -contig {input.contigs} \
            -out {params.output_path}/bin \
            -abund_list {input.depth_file} \
            -min_contig_length {params.minimum_contig_length} \
            -markerset 40 \
            -thread {threads} #-verbose 
        
        mkdir -p {params.output_path}/bins
        
        if grep -q 'Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.' "{params.output_path}/bin.log" \
        || grep -q 'This suggests that the dataset cannot be binned (likely too few and/or small contigs), rather then it actually being an error.' "{params.output_path}/bin.log" \
        || grep -q 'Yielded 0 bins for contig (scaffold) file' "{params.output_path}/bin.log";
        then
            echo "[WARNING]  - Looks like MaxBin2-40 didnt found any prokaryotic bins, this is not a problem and expected for some samples."
            touch "{output.s2b}"
            exit 0
        fi
        
        for FP in {params.output_path}/bin.*.fasta;
        do
            GENOME_SIZE=$(cat $FP | grep -v "^>" | tr -d "\\n" | wc -m)
            echo "[GENOME SIZE] ${{FP}} = ${{GENOME_SIZE}}"
            
            if (( {params.minimum_genome_length} > ${{GENOME_SIZE}} ));
            then
                echo "[COPYING] ${{FP}}"
                ID_GENOME=$(basename ${{FP}} .fasta)
                mv $FP {params.output_path}/bins/${{ID_GENOME}}.fa
            fi
        done
        
        {params.workflow_folder}/scripts/veba/scaffolds_to_bins.py \
            -x fa \
            -i {params.output_path}/bins \
            --bin_prefix {params.output_prefix} \
            > {output.s2b}
        ) &> {log}
        """


checkpoint binning_prokaryotic_dastool:
    input:
        metabat_s2b=rules.binning_prokaryotic_metabat.output.s2b,
        maxbin107_s2b=rules.binning_prokaryotic_maxbin_107.output.s2b,
        maxbin40_s2b=rules.binning_prokaryotic_maxbin_40.output.s2b,
        contigs=get_assembly,
    output:
        bins=directory("samples/{sample}/binning/veba/1_prokaryotic/4_dastool/__DASTool_bins"),
    params:
        workflow_folder=f"{workflow_folder}",
        output_path="samples/{sample}/binning/veba/1_prokaryotic/4_dastool",
        labels="{sample}__METABAT,{sample}__MAXBIN2-107,{sample}__MAXBIN2-40",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/4_dastool.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/4_dastool.log",
    conda:
        "../envs/dastool.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export R_LIBS="$CONDA_PREFIX/lib/R/library"
        
        S2B=$({params.workflow_folder}/scripts/veba/check_scaffolds_to_bins.py \
            -i {input.metabat_s2b},{input.maxbin107_s2b},{input.maxbin40_s2b} \
            -n {params.labels} \
        )
        IFS=" " read -r -a S2B_ARRAY <<< "$S2B"
        
        if [ "${{#S2B_ARRAY[@]}}" -eq 0 ];
        then
            echo "[WARNING] No bins found to combine with DAS_Tool. Skipping."
            mkdir -p "{output.bins}"
            exit 0
        fi
        
        if [ -e "{params.output_path}/manual_skip" ];
        then
            echo "[WARNING] DAS_Tool was skipped manually - Likely becuase of an unrecoverable error in a bad sample."
            mkdir -p "{output.bins}"
            exit 0
        fi
        
        DAS_Tool \
            --bins ${{S2B_ARRAY[0]}} \
            --contigs {input.contigs} \
            --outputbasename {params.output_path}/_ \
            --labels ${{S2B_ARRAY[1]}} \
            --search_engine diamond \
            --score_threshold 0.1 \
            --write_bins \
            --threads {threads} #--debug
        
        if [ ! -d "{output.bins}" ]; then
            echo "[WARNING] {output.bins} does not exist."
            if grep -q 'No bins with bin-score' "{params.output_path}/__DASTool.log";
            then
                echo "[WARNING]  - Looks like DAS_Tool didnt found any prokaryotic bins, this is not a problem and expected for some samples."
                mkdir -p "{output.bins}"
            else
                echo "[ERROR]  - DAS_Tool failed for some reason but didnt return an error. Please check log file."
            fi
        else
            echo "[ERROR]  - DAS_Tool failed to produce bins for some reason but didnt return an error. Please check log file."
        fi
        
        ) &> {log}
        """


rule binning_prokaryotic_whokaryote:
    input:
        fasta=os.path.join(rules.binning_prokaryotic_dastool.output.bins, '{genome}.fa'),
    output:
        fasta="samples/{sample}/binning/veba/1_prokaryotic/5_whokaryote/{genome}/prokaryotes.fasta",
    params:
        output_path="samples/{sample}/binning/veba/1_prokaryotic/5_whokaryote/{genome}",
        minsize=config["veba_prokaryotic"]["minimum_contig_length"],
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/5_whokaryote/{genome}.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/5_whokaryote/{genome}.log",
    conda:
        "../envs/whokaryote.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        whokaryote.py \
            --contigs {input.fasta} \
            --outdir {params.output_path} \
            --minsize {params.minsize} \
            --f \
            --threads {threads} \
        && touch {output.fasta}
        ) &> {log}
        """ # Need to touch output file on sucess since it is not created if we have no prok contigs identified (i.e., is a euk bin)


rule binning_prokaryotic_mdmcleaner:
    input:
        fasta=rules.binning_prokaryotic_whokaryote.output.fasta,
        dbdir=rules.mdmcleaner_download_db.output.dbdir,
    output:
        fasta="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}.cleaned.fa",
    params:
        raw_fasta="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}.fa",
        mdmcleaner_config="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}.mdmcleaner.config",
        output_path="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner",
        output_filtered="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}/{genome}_filtered_kept_contigs.fasta.gz",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner/{genome}.log",
    #conda:
    #    "../envs/mdmcleaner.yaml"
    container:
        "docker://timothystephens/mdmcleaner:0.8.7-TGSv3",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if [ ! -s {input.fasta} ];
        then
            echo '[WARNING] Output from Whokaryote is empty, suggesting that this bin is not prokaryotic. Skipping MDMcleaner.'
            touch {output.fasta}
            exit 0
        fi
        
        cp {input.fasta} {params.raw_fasta}
        echo -e "db_type\\tgtdb\\ndb_basedir\\t{input.dbdir}" > {params.mdmcleaner_config}
        mdmcleaner clean \
            -c {params.mdmcleaner_config} \
            -i {params.raw_fasta} \
            -o {params.output_path} \
            --threads {threads}
        gunzip -c {params.output_filtered} > {output.fasta}
        ) &> {log}
        """


def get_cleaned_prokaryotic_bins(wildcards):
    import os, glob
    
    # [ [sample, genome], [sample, genome], ... ]
    f = glob.glob(
        os.path.join(
            expand(rules.binning_prokaryotic_dastool.output.bins, sample=wildcards.sample)[0],
            "*.fa"
        )
    )
    f = [ [x.split(os.path.sep)[1], os.path.basename(x)] for x in f ]
    f = [ [x[0], x[1].rstrip(".fa")] for x in f ]
    
    outfiles = []
    for x in f:
        outfiles.append(
            expand(rules.binning_prokaryotic_mdmcleaner.output.fasta,
                sample=x[0], genome=x[1]
            )[0]
        )
    
    return outfiles


rule binning_prokaryotic_checkm2:
    input:
        bins=get_cleaned_prokaryotic_bins,
        contigs=get_assembly,
        dbdir=rules.checkm2_download_db.output.dbdir,
        raw_bins=rules.binning_prokaryotic_dastool.output.bins,
    output:
        bins=directory("samples/{sample}/binning/veba/1_prokaryotic/7_checkm2/filtered/genomes"),
        unbinned="samples/{sample}/binning/veba/1_prokaryotic/7_checkm2/filtered/unbinned.fasta",
        quality="samples/{sample}/binning/veba/1_prokaryotic/7_checkm2/filtered/checkm2_results.filtered.tsv",
    params:
        bin_dirs="samples/{sample}/binning/veba/1_prokaryotic/6_mdmcleaner",
        workflow_folder=f"{workflow_folder}",
        output_path="samples/{sample}/binning/veba/1_prokaryotic/7_checkm2",
        tmpdir=temp("tmp/checkm2/samples/{sample}"),
        completeness=config["veba_prokaryotic"]["checkm2_completeness"],
        contamination=config["veba_prokaryotic"]["checkm2_contamination"],
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/7_checkm2.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/7_checkm2.log",
    conda:
        "../envs/checkm2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        rm -fr {params.tmpdir}
        mkdir -p {params.tmpdir}/bins
        
        if ls {params.bin_dirs}/*.cleaned.fa 1> /dev/null 2>&1;
        then
            valid_files=0
            for FA in {params.bin_dirs}/*.cleaned.fa;
            do
                if [ -s "$FA" ];
                then
                    cp "$FA" {params.tmpdir}/bins/$(basename ${{FA%*.cleaned.fa}}).fa
                    valid_files=$((valid_files+1))
                fi
            done
            if [ $valid_files -eq 0 ];
            then
                echo "[WARNING] Only EMPTY bins found. Skipping CheckM2."
                cp {input.contigs} {output.unbinned}
                touch {output.quality}
                mkdir -p {output.bins}
                exit 0
            fi
            
            checkm2 predict \
                -i {params.tmpdir}/bins \
                -o {params.output_path} \
                -t {threads} \
                --force \
                -x fa \
                --tmpdir {params.tmpdir} \
                --database_path {input.dbdir}/CheckM2_database/uniref100.KO.1.dmnd
                
            {params.workflow_folder}/scripts/veba/filter_checkm2_results.py \
                -i {params.output_path}/quality_report.tsv \
                -b {params.tmpdir}/bins \
                -o {params.output_path}/filtered \
                -f {input.contigs} \
                --unbinned \
                --completeness {params.completeness} \
                --contamination {params.contamination} \
                -x fa
        else
            echo "[WARNING] No bins found. Skipping CheckM2."
            cp {input.contigs} {output.unbinned}
            touch {output.quality}
            mkdir -p {output.bins}
        
        fi
        
        ) &> {log}
        """


rule binning_prokaryotic_genome_stats:
    input:
        bins=rules.binning_prokaryotic_checkm2.output.bins,
    output:
        stats="samples/{sample}/binning/veba/1_prokaryotic/7_checkm2/filtered/genome_statistics.tsv",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/1_prokaryotic/8_stats.txt"
    log:
        "logs/samples/{sample}/binning/veba/1_prokaryotic/8_stats.log",
    conda:
        "../envs/seqkit.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if ls {input.bins}/*.fa 1> /dev/null 2>&1;
        then
            seqkit stats \
                -a -b -T -j {threads} \
                {input.bins}/*.fa \
            | python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\\t")' \
            > {output.stats}
        
        else
            echo "[WARNING] No bins found."
            touch {output.stats}
        fi
        ) &> {log}
        """






###############################
####                       ####
####  Binning Eukaryotes   ####
####                       ####
###############################

checkpoint binning_eukaryotic_metabat:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=rules.binning_prokaryotic_checkm2.output.unbinned,
    output:
        s2b="samples/{sample}/binning/veba/2_eukaryotic/1_metabat2/scaffolds_to_bins.tsv",
        bins=directory("samples/{sample}/binning/veba/2_eukaryotic/1_metabat2/bins"),
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_eukaryotic"]["minimum_contig_length"],
        minimum_genome_length=config["veba_eukaryotic"]["minimum_genome_length"],
        output_path="samples/{sample}/binning/veba/2_eukaryotic/1_metabat2",
        output_prefix="{sample}__METABAT2",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/2_eukaryotic/1_metabat2.txt"
    log:
        "logs/samples/{sample}/binning/veba/2_eukaryotic/1_metabat2.log",
    conda:
        "../envs/metabat2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        metabat2 \
            -i {input.contigs} \
            -o {params.output_path}/bins/{params.output_prefix}bin \
            -a {input.depth_file} \
            -m {params.minimum_contig_length} \
            --minClsSize {params.minimum_genome_length} \
            -t {threads} \
            --seed 1 \
            --verbose \
        
        {params.workflow_folder}/scripts/veba/scaffolds_to_bins.py \
            -x fa \
            -i {params.output_path}/bins \
            --bin_prefix {params.output_prefix} \
            > {output.s2b}
        ) &> {log}
        """


rule binning_eukaryotic_whokaryote:
    input:
        fasta=os.path.join(rules.binning_eukaryotic_metabat.output.bins, '{genome}.fa'),
    output:
        fasta="samples/{sample}/binning/veba/2_eukaryotic/2_whokaryote/{genome}/eukaryotes.fasta",
    params:
        output_path="samples/{sample}/binning/veba/2_eukaryotic/2_whokaryote/{genome}",
        minsize=config["veba_eukaryotic"]["minimum_contig_length"],
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/2_eukaryotic/2_whokaryote/{genome}.txt"
    log:
        "logs/samples/{sample}/binning/veba/2_eukaryotic/2_whokaryote/{genome}.log",
    conda:
        "../envs/whokaryote.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        whokaryote.py \
            --contigs {input.fasta} \
            --outdir {params.output_path} \
            --minsize {params.minsize} \
            --f \
            --threads {threads} \
        && touch {output.fasta}
        ) &> {log}
        """ # Need to touch output file on sucess since it is not created if we have no euk contigs identified (i.e., is a prok bin)


rule binning_eukaryotic_busco:
    input:
        fasta=rules.binning_eukaryotic_whokaryote.output.fasta,
        dbdir=rules.busco_download_db.output.dbdir,
    output:
        results=directory("samples/{sample}/binning/veba/2_eukaryotic/3_busco/{genome}"),
    params:
        workflow_folder=f"{workflow_folder}",
        tmp_fasta="samples/{sample}/binning/veba/2_eukaryotic/3_busco/{genome}.fa",
        sample="{sample}",
        completeness=config["veba_eukaryotic"]["busco_completeness"],
        contamination=config["veba_eukaryotic"]["busco_contamination"],
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/2_eukaryotic/3_busco/{genome}.txt"
    log:
        "logs/samples/{sample}/binning/veba/2_eukaryotic/3_busco/{genome}.log",
    #conda:
    #    "../envs/busco.yaml"
    container:
        "docker://timothystephens/busco:6.0.0-TGSv1",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if [ ! -s {input.fasta} ];
        then
            echo '[WARNING] Output from Whokaryote is empty, suggesting that this bin is not eukaryotic. Skipping BUSCO.'
            mkdir -p {output.results}
            exit 0
        fi
        
        export PATH="$CONDA_PREFIX/bin:$PATH"
        export PYTHONPATH="$CONDA_PREFIX/lib/python3.9/site-packages"
        
        mkdir -p {output.results}
        cp {input.fasta} {params.tmp_fasta}
        
        set +eu
        busco \
            --force \
            -i {params.tmp_fasta} \
            -o {output.results} \
            -m genome \
            --auto-lineage-euk \
            -c {threads} \
            --evalue 0.001 \
            --download_path {input.dbdir}
        EXITSTATUS=$?
        
        if [[ "$EXITSTATUS" > 0 ]];
        then
            if grep -q 'EmptyResultsError' "{output.results}/logs/busco.log";
            then
                echo "Looks like BUSCO failed because it found no marker genes, this is not a problem and expected for some bins."
            else
                echo "BUSCO failed. Please check log file."
                exit $EXITSTATUS
            fi
        else
            echo "BUSCO succeeded."
        fi
        
        rm -fr "{output.results}/auto_lineage"
        rm -fr "{output.results}/run_eukaryota_odb10"
        ) &> {log}
        """ # Use custom BUSCO script which handles a BUSCO Exception caused by empty results files nicely (expected sometimes)


def get_cleaned_eukaryotic_bins(wildcards):
    import os, glob
    
    # [ [sample, genome], [sample, genome], ... ]
    f = glob.glob(
        os.path.join(
            expand(rules.binning_eukaryotic_metabat.output.bins, sample=wildcards.sample)[0], 
            "*.fa"
        )
    )
    f = [ [x.split(os.path.sep)[1], os.path.basename(x)] for x in f ]
    f = [ [x[0], x[1].rstrip(".fa")] for x in f ]
    
    o = []
    for x in f:
        o.append(
            expand(rules.binning_eukaryotic_busco.output.results,
                sample=x[0], genome=x[1]
            )[0]
        )
    
    return o


rule binning_eukaryotic_filter:
    input:
        busco=get_cleaned_eukaryotic_bins,
        contigs=rules.binning_prokaryotic_checkm2.output.unbinned,
        raw_bins=rules.binning_eukaryotic_metabat.output.bins,
    output:
        outdir=directory("samples/{sample}/binning/veba/2_eukaryotic/4_filtered"),
        bins=directory("samples/{sample}/binning/veba/2_eukaryotic/4_filtered/genomes"),
        quality="samples/{sample}/binning/veba/2_eukaryotic/4_filtered/busco_results.filtered.tsv",
        tsv="samples/{sample}/binning/veba/2_eukaryotic/4_filtered/busco_results.tsv",
        unbinned="samples/{sample}/binning/veba/2_eukaryotic/4_filtered/unbinned.fasta",
    params:
        workflow_folder=f"{workflow_folder}",
        bins_dir="samples/{sample}/binning/veba/2_eukaryotic/3_busco",
        sample="{sample}",
        minimum_contig_length=config["veba_eukaryotic"]["minimum_contig_length"],
        completeness=config["veba_eukaryotic"]["busco_completeness"],
        contamination=config["veba_eukaryotic"]["busco_contamination"],
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/2_eukaryotic/4_filter.txt"
    log:
        "logs/samples/{sample}/binning/veba/2_eukaryotic/4_filter.log",
    conda:
        "../envs/python.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if [ -z "{input.busco}" ];
        then
            echo "[WARNING] No BUSCO results found. Skipping merging."
            mkdir -p "{output.outdir}" "{output.bins}"
            touch "{output.quality}" "{output.tsv}"
            cp "{input.contigs}" "{output.unbinned}"
            exit 0
        fi
        
        {params.workflow_folder}/scripts/veba/merge_busco_json.py \
            -i {params.bins_dir} \
            -o {output.tsv}
        
        {params.workflow_folder}/scripts/veba/filter_busco_results.py \
            -i {output.tsv} \
            -g {params.bins_dir} \
            -o {output.outdir} \
            -f {input.contigs} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --unbinned
        ) &> {log}
        """


rule binning_eukaryotic_genome_stats:
    input:
        bins=rules.binning_eukaryotic_filter.output.bins,
    output:
        stats="samples/{sample}/binning/veba/2_eukaryotic/4_filtered/genome_statistics.tsv",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/2_eukaryotic/5_stats.txt"
    log:
        "logs/samples/{sample}/binning/veba/2_eukaryotic/5_stats.log",
    conda:
        "../envs/seqkit.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "localrule", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "localrule", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "localrule", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "localrule", "account"),
    shell:
        """
        (
        if ls {input.bins}/*.fa 1> /dev/null 2>&1;
        then
            seqkit stats \
                -a -b -T -j {threads} \
                {input.bins}/*.fa \
              | python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\\t")' \
              > {output.stats}
        else
            echo "[WARNING] No bins found."
            touch {output.stats}
        fi
        ) &> {log}
        """





###############################
####                       ####
####    Binning Viruses    ####
####                       ####
###############################

rule binning_viral_metabat:
    input:
        depth_file=rules.get_metabat_depth_file_combine.output,
        contigs=rules.binning_eukaryotic_filter.output.unbinned,
    output:
        s2b="samples/{sample}/binning/veba/3_viral/1_metabat2/scaffolds_to_bins.tsv",
        merged_bins="samples/{sample}/binning/veba/3_viral/1_metabat2/merged_bins.tsv",
    params:
        workflow_folder=f"{workflow_folder}",
        minimum_contig_length=config["veba_viral"]["minimum_contig_length"],
        minimum_genome_length=config["veba_viral"]["minimum_genome_length"],
        output_path="samples/{sample}/binning/veba/3_viral/1_metabat2",
        output_prefix="{sample}__METABAT2",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/3_viral/1_metabat2.txt"
    log:
        "logs/samples/{sample}/binning/veba/3_viral/1_metabat2.log",
    conda:
        "../envs/metabat2.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        metabat2 \
            -i {input.contigs} \
            -o {params.output_path}/bins/bin \
            -a {input.depth_file} \
            -m {params.minimum_contig_length} \
            --minClsSize {params.minimum_genome_length} \
            -t {threads} \
            --seed 1 \
            --verbose \
        
        {params.workflow_folder}/scripts/veba/scaffolds_to_bins.py \
            -x fa \
            -i {params.output_path}/bins \
            --bin_prefix {params.output_prefix} \
            > {output.s2b}
        
        for MAG in `find {params.output_path}/bins -name "*.fa"`;
        do
            P=$(basename ${{MAG%*.fa}})
            echo ">{params.output_prefix}$P"
            grep -v ">" "$MAG"
        done \
            | seqkit seq -w 0 \
            > {output.merged_bins}
        ) &> {log}
        """


rule binning_viral_genomad:
    input:
        fasta=rules.binning_viral_metabat.output.merged_bins,
        dbdir=rules.genomad_download_db.output.dbdir,
    output:
        virus_summary="samples/{sample}/binning/veba/3_viral/2_genomad/merged_bins_summary/merged_bins_virus_summary.tsv",
        virus_taxonomy="samples/{sample}/binning/veba/3_viral/2_genomad/merged_bins_annotate/merged_bins_taxonomy.tsv",
        plasmid_summary="samples/{sample}/binning/veba/3_viral/2_genomad/merged_bins_summary/merged_bins_plasmid_summary.tsv",
    params:
        results="samples/{sample}/binning/veba/3_viral/2_genomad",
        workflow_folder=f"{workflow_folder}",
        sample="{sample}",
        empty_virus_summary=f"{workflow_folder}/data/virus_summary.tsv",
        empty_virus_taxonomy=f"{workflow_folder}/data/virus_taxonomy.tsv",
        empty_plasmid_summary=f"{workflow_folder}/data/plasmid_summary.tsv",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/3_viral/2_genomad.txt"
    log:
        "logs/samples/{sample}/binning/veba/3_viral/2_genomad.log",
    #conda:
    #    "../envs/genomad.yaml"
    container:
        "docker://antoniopcamargo/genomad:1.11.0",
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if [ -s {input.fasta} ];
        then
          export PATH="/opt/conda/bin:$PATH"
          genomad end-to-end \
              --cleanup \
              --threads {threads} \
              --verbose \
              --enable-score-calibration \
              --disable-find-proviruses \
              --sensitivity 4.0 \
              --splits 0 \
              --composition auto \
              --min-score 0.7 \
              --max-fdr 1.0 \
              --min-plasmid-marker-enrichment -100 \
              --min-virus-marker-enrichment -100 \
              --min-plasmid-hallmarks 0 \
              --min-virus-hallmarks 0 \
              --max-uscg 100 \
              {input.fasta} \
              {params.results} \
              {input.dbdir}/genomad_db
        else
          echo '[WARNING] Output from MetaBAT2 is empty, suggesting that no bins were present in the remaining scaffolds. Skipping geNomad and creating empty output files for downstream analysis.'
          mkdir -p "{params.sample}/binning/veba/3_viral/2_genomad/merged_bins_summary"
          mkdir -p "{params.sample}/binning/veba/3_viral/2_genomad/merged_bins_annotate"
          cat "{params.empty_virus_summary}"   > "{output.virus_summary}"
          cat "{params.empty_virus_taxonomy}"  > "{output.virus_taxonomy}"
          cat "{params.empty_plasmid_summary}" > "{output.plasmid_summary}"
        fi
        ) &> {log}
        """


rule binning_viral_filter:
    input:
        contigs=rules.binning_eukaryotic_filter.output.unbinned,
        s2b=rules.binning_viral_metabat.output.s2b,
        virus_summary=rules.binning_viral_genomad.output.virus_summary,
        virus_taxonomy=rules.binning_viral_genomad.output.virus_taxonomy,
        plasmid_summary=rules.binning_viral_genomad.output.plasmid_summary,
    output:
        plasmid_bins=directory("samples/{sample}/binning/veba/3_viral/2_genomad/filtered_plasmid_bins/genomes"),
        plasmid_completness="samples/{sample}/binning/veba/3_viral/2_genomad/filtered_plasmid_bins/genomad_results.filtered.tsv",
        viral_bins=directory("samples/{sample}/binning/veba/3_viral/2_genomad/filtered_viral_bins/genomes"),
        viral_completness="samples/{sample}/binning/veba/3_viral/2_genomad/filtered_viral_bins/genomad_results.filtered.tsv",
        unbinned="samples/{sample}/binning/veba/3_viral/2_genomad/unbinned.fasta",
    params:
        plasmid_output=directory("samples/{sample}/binning/veba/3_viral/2_genomad/filtered_plasmid_bins"),
        viral_output=directory("samples/{sample}/binning/veba/3_viral/2_genomad/filtered_viral_bins"),
        workflow_folder=f"{workflow_folder}",
        outdir="samples/{sample}/binning/veba/3_viral/2_genomad",
        sample="{sample}",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/3_viral/3_filter.txt"
    log:
        "logs/samples/{sample}/binning/veba/3_viral/3_filter.log",
    conda:
        "../envs/virus_filter.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        echo "# Filtering Virus Results"
        {params.workflow_folder}/scripts/veba/filter_genomad_results.py \
            --fasta {input.contigs} \
            --scaffolds_to_bins {input.s2b} \
            --genomad_results {input.virus_summary} \
            --genomad_virus_taxonomy {input.virus_taxonomy} \
            --prefix {params.sample}__VIRUS__ \
            --output_directory {params.viral_output}
        
        echo "# Filtering Plasmid Results"
        {params.workflow_folder}/scripts/veba/filter_genomad_results.py \
            --fasta {input.contigs} \
            --scaffolds_to_bins {input.s2b} \
            --genomad_results {input.plasmid_summary} \
            --prefix {params.sample}__PLASMID__ \
            --output_directory {params.plasmid_output}
        
        echo "# Getting unbinned sequences"
        cat {params.viral_output}/unbinned.list {params.plasmid_output}/unbinned.list \
            | sort | uniq -d \
            > {params.outdir}/unbinned.list \
          && \
        cat {input.contigs} \
            | seqkit grep \
                --pattern-file {params.outdir}/unbinned.list \
            > {output.unbinned}
        ) &> {log}
        """


rule binning_viral_genome_stats:
    input:
        plasmid_bins=rules.binning_viral_filter.output.plasmid_bins,
        viral_bins=rules.binning_viral_filter.output.viral_bins,
    output:
        plasmid_stats="samples/{sample}/binning/veba/3_viral/2_genomad/filtered_plasmid_bins/genome_statistics.tsv",
        viral_stats="samples/{sample}/binning/veba/3_viral/2_genomad/filtered_viral_bins/genome_statistics.tsv",
    benchmark:
        "benchmarks/samples/{sample}/binning/veba/3_viral/4_stats.txt"
    log:
        "logs/samples/{sample}/binning/veba/3_viral/4_stats.log",
    conda:
        "../envs/seqkit.yaml"
    threads: lambda wc: get_resource(wc, None, 1, "binning", "threads")
    resources:
        mem_mb=lambda wc, input, attempt: get_resource(wc, input, attempt, "binning", "mem_mb"),
        partition=lambda wildcards, input, attempt: get_queue(input, attempt, "binning", "partition"),
        account=lambda wildcards, input, attempt: get_resource(wildcards, input, attempt, "binning", "account"),
    shell:
        """
        (
        if ls {input.plasmid_bins}/*.fa 1> /dev/null 2>&1;
        then
            seqkit stats \
                -a -b -T -j {threads} \
                {input.plasmid_bins}/*.fa \
              | python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\\t")' \
              > {output.plasmid_stats}
        else
            echo "[WARNING] No Plasmid bins found."
            touch {output.plasmid_stats}
        fi
        
        if ls {input.viral_bins}/*.fa 1> /dev/null 2>&1;
        then
            seqkit stats \
                -a -b -T -j {threads} \
                {input.viral_bins}/*.fa \
              | python -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep="\t", index_col=0); df.index = df.index.map(lambda x: x[:-3]); df.to_csv(sys.stdout, sep="\\t")' \
              > {output.viral_stats}
        else
            echo "[WARNING] No Viral bins found."
            touch {output.viral_stats}
        fi
        ) &> {log}
        """
