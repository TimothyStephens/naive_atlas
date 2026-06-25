import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def gff_to_prodigal(gff_files, fna_files, out_file, skipped_out_file, map_out_file):
    # 1. Load all contig sequences from all FNA files into one master dictionary
    contigs = {}
    for fna_file in fna_files:
        print(f"Loading FASTA: {fna_file}")
        for record in SeqIO.parse(fna_file, "fasta"):
            if record.id in contigs:
                print(f"WARNING: Duplicate contig ID '{record.id}' found. Overwriting previous entry.")
            contigs[record.id] = record

    # 2. Parse all GFF files and group features
    features = defaultdict(list)
    
    # Exclude structural/parent features that are redundant to the actual sequence blocks
    redundant_features = {"gene", "mRNA", "exon", "CRISPR"}
    
    for gff_file in gff_files:
        print(f"Parsing GFF: {gff_file}")
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                    
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                    
                scaff = parts[0]
                feature_type = parts[2]
                
                # Skip the redundant features
                if feature_type in redundant_features:
                    continue
                    
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes_str = parts[8]
                
                # Extract attributes into a dictionary
                attr_dict = {}
                for attr in attributes_str.split(";"):
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attr_dict[k] = v
                
                # GROUPING LOGIC:
                # CDS uses Parent to catch multi-exons. 
                # Non-CDS uses its own ID so distinct features aren't falsely merged.
                if feature_type == "CDS":
                    group_id = attr_dict.get("Parent", attr_dict.get("gene_id", attr_dict.get("ID", "unknown")))
                else:
                    group_id = attr_dict.get("ID", attr_dict.get("Parent", f"{feature_type}_{scaff}_{start}_{end}"))
                    
                features[group_id].append((scaff, start, end, strand, feature_type))
            
    # 3. Process features, check for multi-segments, format as Prodigal, and map IDs
    records_to_write = []
    skipped_features = []
    feature_mapping = []
    scaff_counters = defaultdict(int)
    
    for group_id, loc_list in features.items():
        # Multi-segment check: Skip and log to the skipped warning list
        if len(loc_list) > 1:
            skipped_features.append(group_id)
            continue
        
        scaff, start, end, strand, feature_type = loc_list[0]
        
        if scaff not in contigs:
            raise ValueError(f"ERROR: Scaffold '{scaff}' found in GFF for '{group_id}' but not in any FNA files.")
            
        # Keep a running counter to generate Prodigal IDs
        scaff_counters[scaff] += 1
        
        # APPEND FEATURE TYPE DIRECTLY TO NUMBER (no underscore)
        # E.g. scaff_1CDS. This prevents external tools from splitting the scaffold name incorrectly.
        prodigal_id = f"{scaff}_{scaff_counters[scaff]}{feature_type}"
        
        # Store the mapping 
        feature_mapping.append((group_id, prodigal_id, feature_type))
        
        # Convert strand notation (+/- to 1/-1). Default to 1 for unstranded features.
        strand_num = "-1" if strand == "-" else "1"
        
        # GFF coordinates are 1-based, inclusive. Python slicing is 0-based, exclusive.
        seq_slice = contigs[scaff].seq[start - 1 : end]
        
        if strand == "-":
            seq_slice = seq_slice.reverse_complement()
            
        # Sanity check against the specific length validation in your downstream tool
        if (end - start) + 1 != len(seq_slice):
            raise ValueError(
                f"ERROR: Length mismatch for {prodigal_id}; "
                f"start={start} end={end} length={len(seq_slice)}"
            )
            
        # Format the description exactly as parse_prodigal_genes expects it.
        desc = f"# {start} # {end} # {strand_num} # ID={prodigal_id};partial=00;type={feature_type}"
        
        rec = SeqRecord(
            seq_slice,
            id=prodigal_id,
            description=desc
        )
        records_to_write.append(rec)
        
    # 4. Write out the Prodigal formatted FASTA
    print(f"Writing {len(records_to_write)} features to {out_file}")
    SeqIO.write(records_to_write, out_file, "fasta")

    # 5. Write out the mapping file (TSV format)
    print(f"Writing feature ID mapping to {map_out_file}")
    with open(map_out_file, 'w') as f:
        f.write("original_gff_id\tprodigal_formatted_id\tfeature_type\n")
        for orig_id, new_id, f_type in feature_mapping:
            f.write(f"{orig_id}\t{new_id}\t{f_type}\n")

    # 6. Write out the skipped features warning file
    print(f"Skipped {len(skipped_features)} multi-segment features. Writing list to {skipped_out_file}")
    with open(skipped_out_file, 'w') as f:
        f.write("# List of multi-segment features skipped during Prodigal format conversion\n")
        for skipped in skipped_features:
            f.write(f"{skipped}\n")


# --- Snakemake Integration ---
if __name__ == "__main__":
    try:
        gff_in = snakemake.params.gff if isinstance(snakemake.params.gff, list) else [snakemake.params.gff]
        fna_in = snakemake.params.fna if isinstance(snakemake.params.fna, list) else [snakemake.params.fna]
        
        fna_out = snakemake.output.fna
        skipped_out = snakemake.output.skipped
        map_out = snakemake.output.mapping
        
        gff_to_prodigal(gff_in, fna_in, fna_out, skipped_out, map_out)
        
    except NameError:
        print("This script is designed to be run within a Snakemake rule.")
        sys.exit(1)
