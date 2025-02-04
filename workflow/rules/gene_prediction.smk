


#################################
####                         ####
####      Predict Genes      ####
####                         ####
#################################

def get_genomes_for_gene_prediction(lineage):
    genome_dir = 'genomes/genomes'
    
    if lineage == 'prokaryotic':
        fasta_files = glob(os.path.join(genome_dir, "MAG_prokaryotic_*.fa"))
        if len(fasta_files) == 0:
            print(f"No Prokaryotic genomes found with fa extension in {genome_dir} ")
            return([])
        
        genomes = []
        for file_path in fasta_files:
            file_name = os.path.basename(file_path)  # "file.txt"
            file_name_without_ext = os.path.splitext(file_name)[0]  # "file"
            genomes.append(file_name_without_ext)
        return(genomes)
    
    if lineage == 'eukaryotic':
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
    
    if lineage == 'viral':
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
        fna="Predict_Genes/genomes/bacteria/{genome}.fna",
        gff="Predict_Genes/genomes/bacteria/{genome}.gff",
        rrna="Predict_Genes/genomes/bacteria/{genome}.rRNA.fna",
        trna="Predict_Genes/genomes/bacteria/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
        wd="Predict_Genes/genomes/bacteria",
        genome="{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/bacteria/{genome}.txt"
    log:
        "logs/gene_prediction/bacteria/{genome}.txt",
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
            --output {params.wd} --force \
            --meta \
            --threads {threads} \
            --keep-contig-headers \
            {input.fasta}
        ) &> {log}
        """


rule gene_prediction_archaea:
    input:
        fasta="genomes/genomes/{genome}.fa",
    output:
        faa="Predict_Genes/genomes/archaea/{genome}.faa",
        fna="Predict_Genes/genomes/archaea/{genome}.fna",
        gff="Predict_Genes/genomes/archaea/{genome}.gff",
        rrna="Predict_Genes/genomes/archaea/{genome}.rRNA.fna",
        trna="Predict_Genes/genomes/archaea/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
    benchmark:
        "logs/benchmarks/gene_prediction/archaea/{genome}.txt"
    log:
        "logs/gene_prediction/archaea/{genome}.txt",
    conda:
        "../envs/gene_prediction_archaea.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        prokka --outdir tmp --force --prefix MAG_prokaryotic_03 --cpus 8 --addgenes --addmrna --kingdom Archaea  MAG_prokaryotic_03.fa
        ) &> {log}
        """


rule gene_prediction_virus:
    input:
        fasta="genomes/genomes/{genome}.fa",
    output:
        faa="Predict_Genes/genomes/virus/{genome}.faa",
        fna="Predict_Genes/genomes/virus/{genome}.fna",
        gff="Predict_Genes/genomes/virus/{genome}.gff",
        rrna="Predict_Genes/genomes/virus/{genome}.rRNA.fna",
        trna="Predict_Genes/genomes/virus/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
    benchmark:
        "logs/benchmarks/gene_prediction/virus/{genome}.txt"
    log:
        "logs/gene_prediction/virus/{genome}.txt",
    conda:
        "../envs/gene_prediction_virus.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        prokka --outdir tmp --force --prefix MAG_viral_001      --cpus 8 --addgenes --addmrna --kingdom Viruses  MAG_viral_001.fa
        ) &> {log}
        """


rule gene_prediction_plasmid:
    input:
        fasta="genomes/genomes/{genome}.fa",
    output:
        faa="Predict_Genes/genomes/plasmid/{genome}.faa",
        fna="Predict_Genes/genomes/plasmid/{genome}.fna",
        gff="Predict_Genes/genomes/plasmid/{genome}.gff",
        rrna="Predict_Genes/genomes/plasmid/{genome}.rRNA.fna",
        trna="Predict_Genes/genomes/plasmid/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
    benchmark:
        "logs/benchmarks/gene_prediction/plasmid/{genome}.txt"
    log:
        "logs/gene_prediction/plasmid/{genome}.txt",
    conda:
        "../envs/gene_prediction_plasmid.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_memory"],
        time=config["simplejob_runtime"],
    shell:
        """
        (
        bakta --db tmp/db/db --prefix MAG_plasmid_01     --output tmp --force --meta --threads 8 MAG_plasmid_01.fa
        ) &> {log}
        """






rule gene_prediction_eukaryotic:
    input:
        fasta=expand("genomes/genomes/{genome}.fa", genomes=get_genomes_for_gene_prediction('eukaryotic')),
        dbdir=rules.microeukaryotic_mmseqs2_db.output.dbdir,
    output:
        faa="genomes/genomes/{genome}.faa",
        fna="genomes/genomes/{genome}.fna",
        gff="genomes/genomes/{genome}.gff",
        rrna="genomes/genomes/{genome}.rRNA.fna",
        trna="genomes/genomes/{genome}.tRNA.fna",
    params:
        workflow_folder=f"{workflow_folder}",
        outdir="Gene_Prediction/eukaryotic/{genome}",
    benchmark:
        "logs/benchmarks/gene_prediction/{genome}.txt"
    log:
        "logs/gene_prediction/{genome}.txt",
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






