import pandas as pd

import logging

logger = logging.getLogger(__file__)


def load_sample_table(sample_table="samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    validate_bingroup_size(sampleTable)
    return sampleTable


def validate_sample_table(sampleTable):
    Expected_Headers = ["Reads_raw_R1", "Reads_raw_R2", "Reads_raw_Long", "Normalize_reads_before_assembly", "Error_correction_before_assembly", "Assembler", "Bin_group"]
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
            logger.error(f"expect '{h}' to be found in samples.tsv")
            exit(1)
        #elif sampleTable[h].isnull().any():
        #    logger.error(f"Found empty values in the sample table column '{h}'")
        #    exit(1)

    if not sampleTable.index.is_unique:
        duplicated_samples = ", ".join(sampleTable.index.duplicated())
        logger.error(
            f"Expect Samples to be unique. Found {duplicated_samples} more than once"
        )
        exit(1)

    if sampleTable.index.str.match("^\d").any():
        logger.error(
            f"Sample names shouldn't start with a digit. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )
        exit(1)

    if sampleTable.index.str.contains("_").any():
        logger.error(
            f"Sample names shouldn't contain underscores. This can lead to incompatibilities. \n {list(sampleTable.index)}"
        )
        exit(1)

    if sampleTable.index.str.count("-").max() > 1:
        logger.error(
            f"Sample names shouldn't have more than one hypo '-'. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )
        exit(1)

    ### Validate Bin_group

    if sampleTable.Bin_group.isnull().any():
        logger.warning(f"Found empty values in the sample table column 'Bin_group'")

    if sampleTable.Bin_group.str.contains("_").any():
        logger.error(
            f"Bin_group names shouldn't contain underscores. This can lead to incompatibilities. \n {list(sampleTable.Bin_group)}"
        )
        exit(1)

    if sampleTable.Bin_group.str.contains("-").any():
        logger.error(
            f"Bin_group names shouldn't contain hypos '-'. This can lead to incompatibilities.\n {list(sampleTable.Bin_group)}"
        )
        exit(1)
    
    ### Validate Assembler
    sampleTable['Assembler'] = sampleTable['Assembler'].str.lower()
    unknown_assemblers = [x for x in sampleTable.Assembler.unique() if not x in ['megahit', 'spades', 'flye']]
    if unknown_assemblers:
        logger.error(
            f"Assembler found that is not part of allowable options: 'megahit', 'spades', and 'flye'.\nBad assembler options: {unknown_assemblers}"
        )
        exit(1)
    
    ### Enforce Bool
    sampleTable['Normalize_reads_before_assembly']  = sampleTable['Normalize_reads_before_assembly'].map(
            {'True': True, 'T': True, 'true': True, '1': True, 'False': False, 'F': False, 'false': False, '0': False}
        ).fillna(False).astype(bool)
    sampleTable['Error_correction_before_assembly'] = sampleTable['Error_correction_before_assembly'].map(
            {'True': True, 'T': True, 'true': True, '1': True, 'False': False, 'F': False, 'false': False, '0': False}
        ).fillna(False).astype(bool)


class BinGroupSizeError(Exception):
    """
    Exception with Bingroupsize
    """

    def __init__(self, message):
        super(BinGroupSizeError, self).__init__(message)


def validate_bingroup_size(sampleTable):
    bin_group_sizes = sampleTable.Bin_group.value_counts()

    max_bin_group_size = bin_group_sizes.max()

    warn_message = (
        "Co-binning uses cross-mapping which scales quadratically."
        f"You have a bingroup with {max_bin_group_size} samples, which already leads to {max_bin_group_size*max_bin_group_size} cross-mappings."
    )

    if max_bin_group_size > 50:
        logger.error(
            warn_message
            + f"Max bin group size of {max_bin_group_size} is too much . Please split your samples into smaller groups."
        )
        BinGroupSizeError("Bin_group too large")

    if max_bin_group_size > 15:
        logger.warning(
            warn_message
            + f"Max bin group size of {max_bin_group_size} might be too much for cross-mapping. Consider spliting your samples into smaller groups."
        )

    elif max_bin_group_size == 1:
        logger.warning(
            "You have only one sample per bingroup. Will use this information but your bins might be less accurate."
        )


