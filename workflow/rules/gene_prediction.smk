


#################################
####                         ####
####  Predict Genes (Genome) ####
####                         ####
#################################

def get_genomes_for_gene_prediction(lineage):
    import pandas as pd
    
    genome_dir = 'genomes/genomes'
    
    if lineage == 'bacteria':
        fasta_files = glob(os.path.join(genome_dir, "MAG_prokaryotic_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Prokaryotic genomes found with fa extension in {genome_dir} ")
            return([])
        
        file_name = "genomes/annotations/genomes/taxonomy/gtdb_taxonomy.tsv"
        annot = pd.read_table(file_name, sep='\t', index_col=0)
        
        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            if annot.loc[file_name_without_ext, "Domain"] == "Bacteria":
                genomes.append(file_name_without_ext)
        return(genomes)
    
    if lineage == 'archaea':
        fasta_files = glob(os.path.join(genome_dir, "MAG_prokaryotic_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Prokaryotic genomes found with fa extension in {genome_dir} ")
            return([])

        file_name = "genomes/annotations/genomes/taxonomy/gtdb_taxonomy.tsv"
        annot = pd.read_table(file_name, sep='\t', index_col=0)

        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            if annot.loc[file_name_without_ext, "Domain"] == "Archaea":
                genomes.append(file_name_without_ext)
        return(genomes)
    
    if lineage == 'eukaryote':
        fasta_files = glob(os.path.join(genome_dir, "MAG_eukaryotic_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Eukaryotic genomes found with fa extension in {genome_dir} ")
            return([])
        
        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            genomes.append(file_name_without_ext)
        return(genomes)
    
    if lineage == 'virus':
        fasta_files = glob(os.path.join(genome_dir, "MAG_viral_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Viral genomes found with fa extension in {genome_dir} ")
            return([])
        
        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            genomes.append(file_name_without_ext)
        return(genomes)
    
    if lineage == 'plasmid':
        fasta_files = glob(os.path.join(genome_dir, "MAG_plasmid_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Plastid genomes found with fa extension in {genome_dir} ")
            return([])
        
        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            genomes.append(file_name_without_ext)
        return(genomes)


rule gene_prediction_bacteria:
    input:
        fasta="genomes/genomes/{genome}.fa",
        dbdir=rules.bakta_download_db.output.dbdir,
    output:
        faa="Predict_Genes/genomes/bacteria/{genome}.faa",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/genomes/bacteria",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/genomes/bacteria/{genome}.txt"
    log:
        "logs/gene_prediction/genomes/bacteria/{genome}.txt",
    conda:
        "../envs/gene_prediction_bacteria.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        rm -fr {params.wd}/{params.genome}*
        
        bakta \
            --db {input.dbdir} \
            --prefix {params.genome} \
            --locus {params.genome} \
            --locus-tag {params.genome} \
            --output {params.wd} --force \
            --meta \
            --keep-contig-headers \
            --threads {threads} \
            {input.fasta}
        ) &> {log}
        """


rule gene_prediction_archaea:
    input:
        fasta="genomes/genomes/{genome}.fa",
    output:
        faa="Predict_Genes/genomes/archaea/{genome}.faa",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/genomes/archaea",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/genomes/archaea/{genome}.txt"
    log:
        "logs/gene_prediction/genomes/archaea/{genome}.txt",
    conda:
        "../envs/gene_prediction_archaea.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        prokka \
            --outdir {params.wd} --force \
            --prefix {params.genome} \
            --locustag {params.genome} \
            --cpus {threads} \
            --addgenes --addmrna --metagenome \
            --kingdom Archaea \
            {input.fasta}
        ) &> {log}
        """


rule gene_prediction_virus:
    input:
        fasta="genomes/genomes/{genome}.fa",
    output:
        faa="Predict_Genes/genomes/virus/{genome}.faa",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/genomes/virus",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/genomes/virus/{genome}.txt"
    log:
        "logs/gene_prediction/genomes/virus/{genome}.txt",
    conda:
        "../envs/gene_prediction_virus.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        export PERL5LIB="$CONDA_PREFIX/lib/site_perl/5.26.2"
        prokka \
            --outdir {params.wd} --force \
            --prefix {params.genome} \
            --locustag {params.genome} \
            --cpus {threads} \
            --addgenes --addmrna --metagenome \
            --kingdom Viruses \
            {input.fasta}
        ) &> {log}
        """


rule gene_prediction_plasmid:
    input:
        fasta="genomes/genomes/{genome}.fa",
        dbdir=rules.bakta_download_db.output.dbdir,
    output:
        faa="Predict_Genes/genomes/plasmid/{genome}.faa",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/genomes/plasmid",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/genomes/plasmid/{genome}.txt"
    log:
        "logs/gene_prediction/genomes/plasmid/{genome}.txt",
    conda:
        "../envs/gene_prediction_plasmid.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        bakta \
            --db {input.dbdir} \
            --prefix {params.genome} \
            --locus {params.genome} \
            --locus-tag {params.genome} \
            --output {params.wd} --force \
            --meta \
            --keep-contig-headers \
            --threads {threads} \
            {input.fasta}
        ) &> {log}
        """


rule gene_prediction_eukaryote:
    input:
        fasta="genomes/genomes/{genome}.fa",
        dbdir=rules.microeukaryotic_mmseqs2_db.output.dbdir,
    output:
        faa="Predict_Genes/genomes/eukaryotes/{genome}.faa",
        fna="Predict_Genes/genomes/eukaryotes/{genome}.fna",
        gff="Predict_Genes/genomes/eukaryotes/{genome}.gff",
        rrna="Predict_Genes/genomes/eukaryotes/{genome}.rRNA.fna",
        trna="Predict_Genes/genomes/eukaryotes/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
        outdir="Gene_Prediction/eukaryotic/{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/genomes/{genome}.txt"
    log:
        "logs/gene_prediction/genomes/{genome}.txt",
    conda:
        "../envs/gene_prediction_eukaryotic.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        mkdir -p {params.outdir}
        {params.workflow_folder}/scripts/veba/eukaryotic_gene_modeling_wrapper.py \
            --fasta {input.fasta} \
            --metaeuk_database {input.dbdir}/MicroEuk50 \
            --metaeuk_split_memory_limit 36G \
            -o {params.outdir} \
            -p {threads} \
            --metaeuk_sensitivity 4.0 \
            --metaeuk_evalue 0.01 \
            --pyrodigal_minimum_gene_length 90 \
            --pyrodigal_minimum_edge_gene_length 60 \
            --pyrodigal_maximum_gene_overlap_length 60 \
            --pyrodigal_mitochondrial_genetic_code 4 \
            --pyrodigal_plastid_genetic_code 11 \
            --barrnap_length_cutoff 0.8 \
            --barrnap_reject 0.25 \
            --barrnap_evalue 1e-06 \
            --trnascan_mitochondrial_searchmode='-O' \
            --trnascan_plastid_searchmode='-O'
        ) &> {log}
        """


rule move_genome_predicted_genes:
    input:
        bacteria=expand(rules.gene_prediction_bacteria.output.faa,
                    genome=get_genomes_for_gene_prediction("bacteria")),
        archaea=expand(rules.gene_prediction_archaea.output.faa,
                    genome=get_genomes_for_gene_prediction("archaea")),
        eukaryote=expand(rules.gene_prediction_eukaryote.output.faa,
                    genome=get_genomes_for_gene_prediction("eukaryote")),
        viral=expand(rules.gene_prediction_virus.output.faa,
                    genome=get_genomes_for_gene_prediction("virus")),
        plasmid=expand(rules.gene_prediction_plasmid.output.faa,
                    genome=get_genomes_for_gene_prediction("plasmid")),
    output:
        outdir=directory("genomes/genes/genomes"),
    params:
        workflow_folder=f"{workflow_folder}",
        prokaryotic_stats="genomes/genes/genomes/MAG_prokaryotic.gene_stats.tsv",
        eukaryotic_stats="genomes/genes/genomes/MAG_eukaryotic.gene_stats.tsv",
        viral_stats="genomes/genes/genomes/MAG_viral.gene_stats.tsv",
        plasmid_stats="genomes/genes/genomes/MAG_plasmid.gene_stats.tsv",
    log:
        "logs/gene_prediction/genomes/move_predicted_genes.txt",
    conda:
        "../envs/python.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        rm -fr {output.outdir}; mkdir -p {output.outdir}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_PROKKA.py \
            -i {input.archaea} \
            -o {output.outdir} \
            -s {params.prokaryotic_stats}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_BAKTA.py \
            -i {input.bacteria} \
            -o {output.outdir} \
            -s {params.prokaryotic_stats}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_PROKKA.py \
            -i {input.viral} \
            -o {output.outdir} \
            -s {params.viral_stats}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_BAKTA.py \
            -i {input.plasmid} \
            -o {output.outdir} \
            -s {params.plasmid_stats}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_EUK.py \
            -i {input.eukaryote} \
            -o {output.outdir} \
            -s {params.eukaryotic_stats}
        
        ) &> {log}
        """





#################################
####                         ####
#### Predict Genes (Unbinned)####
####                         ####
#################################

rule gene_prediction_unbinned:
    input:
        fasta="genomes/unbinned/{genome}.fa",
    output:
        faa="Predict_Genes/unbinned/{genome}.faa",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/unbinned",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/unbinned/{genome}.txt"
    log:
        "logs/gene_prediction/unbinned/{genome}.txt",
    conda:
        "../envs/prodigal.yaml"
    threads: 1
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        rm -fr {params.wd}/{params.genome}*
        
        prodigal \
            -a {params.wd}/{params.genome}.faa \
            -d {params.wd}/{params.genome}.fna \
            -f gff \
            -o {params.wd}/{params.genome}.gff3 \
            -p meta \
            < {input.fasta}
        ) &> {log}
        """

def get_all_unbinned_genes(wildcards):
    all_genomes = get_all_unbinned(wildcards)
    return(
        expand(rules.gene_prediction_unbinned.output.faa,
               genome=all_genomes)
    )

rule move_unbinned_predicted_genes:
    input:
        unbinned=get_all_unbinned_genes,
    output:
        outdir=directory("genomes/genes/unbinned"),
    params:
        workflow_folder=f"{workflow_folder}",
        unbinned_stats="genomes/genes/unbinned/Unbinned.gene_stats.tsv",
    log:
        "logs/gene_prediction/unbinned/move_predicted_genes.txt",
    conda:
        "../envs/python.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        rm -fr {output.outdir}; mkdir -p {output.outdir}
        
        {params.workflow_folder}/scripts/veba/prepare_predicted_genes_from_PRODIGAL.py \
            -i {input.unbinned} \
            -o {output.outdir} \
            -s {params.unbinned_stats}
        
        ) &> {log}
        """


