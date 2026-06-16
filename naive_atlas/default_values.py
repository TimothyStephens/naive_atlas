
#### Hardcoded paths
spaces_found = 0

# Use data/adapters.fa unless overwritten by the user
ADAPTERS = config.get("preprocess_adapters", os.path.join(workflow_folder, "data/adapters.fa"))
ADAPTERS = os.path.realpath(ADAPTERS)
if " " in ADAPTERS:
    logger.error(f"The adapter path contains spaces: '{ADAPTERS}'. Paths with spaces are not supported.")
    spaces_found += 1

# Use PhiX (unless overwritten by the user) and append user added references
CONTAMINANT_REFERENCES = {}
if config.get("contaminant_references_include_phiX"):
    CONTAMINANT_REFERENCES["PhiX"] = os.path.join(workflow_folder, "data/phiX174_virus.fa")
CONTAMINANT_REFERENCES = CONTAMINANT_REFERENCES | config.get("contaminant_references", {})

for name, path in CONTAMINANT_REFERENCES.items():
    # Handle if path is wrapped in a list/tuple (e.g. metatranscriptome rRNA)
    actual_path = path[0] if isinstance(path, (list, tuple)) else path
    if actual_path:
        real_path = os.path.realpath(actual_path)
        if " " in real_path:
            logger.error(f"Contaminant reference '{name}' path contains spaces: '{real_path}'. Paths with spaces are not supported.")
            spaces_found += 1
        CONTAMINANT_REFERENCES[name] = real_path

if spaces_found > 0:
    sys.exit(1)

logger.debug(f"Adapter file being used: {ADAPTERS}")
logger.debug(f"contaminant reference files being used: {CONTAMINANT_REFERENCES}")



EGGNOG_HEADER = [
    "Query",
    "Seed",
    "Seed_evalue",
    "Seed_Score",
    "eggNOG",
    "max_annot_lvl",
    "COG_cat",
    "Description",
    "Name",
    "GO_terms",
    "EC",
    "KO",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]
