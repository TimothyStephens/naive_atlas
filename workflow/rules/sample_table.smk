from naive_atlas.sample_table import load_sample_table, validate_bingroup_size

sampleTable = load_sample_table()

def io_params_for_tadpole(io, key="in"):
    """This function generates the input flag needed for bbwrap/tadpole for all cases
    possible for get_quality_controlled_reads.

    params:
        io  input or output element from snakemake
        key 'in' or 'out'

        if io contains attributes:
            R1 -> in={R1/LR}
            R1,R2 -> in1={R1} in2={R2}

    """
    N = len(io)
    if N == 1:
        flag = f"{key}1={io[0]}"
    elif N == 2:
        flag = f"{key}1={io[0]} {key}2={io[1]}"

    else:
        logger.error(
            (
                "File input/output expectation is one of: "
                "1 file = single-end/ interleaved paired-end "
                "2 files = R1,R2, or"
                "got: {n} files:\n{}"
            ).format("\n".join(io), n=len(io))
        )
        sys.exit(1)
    return flag


def input_params_for_bbwrap(input):
    if len(input) == 3:
        return f"in1={input[0]},{input[2]} in2={input[1]},null"
    else:
        return io_params_for_tadpole(input)


SAMPLES = sampleTable.index.values


# GROUPS = sampleTable.Bin_group.unique()
def get_alls_samples_of_group(wildcards):
    group_of_sample = sampleTable.loc[wildcards.sample, "Bin_group"]

    return list(sampleTable.loc[sampleTable.Bin_group == group_of_sample].index)


class FileNotInSampleTableException(Exception):
    """
    Exception with sampleTable
    """

    def __init__(self, message):
        super(FileNotInSampleTableException, self).__init__(message)


def get_files_from_sampleTable(sample, Headers):
    """
    Function that gets some filenames form the sampleTable for a given sample and Headers.
    It checks various possibilities for errors and throws either a
    FileNotInSampleTableException or a IOError, when something went really wrong.
    """
    
    if not (sample in sampleTable.index):
        raise IOError(f"Sample name {sample} is not in sampleTable")
    
    if not (sample in sampleTable.index):
        raise IOError(f"Sample name {sample} is not in sampleTable")

    Error_details = f"\nsample: {sample}\nFiles: {Headers}"

    if type(Headers) == str:
        Headers = [Headers]

    NheadersFound = sampleTable.columns.isin(Headers).sum()

    if NheadersFound == 0:
        raise FileNotInSampleTableException(
            f"None of the Files are in sampleTable, they should be added to the sampleTable later in the workflow"
            + Error_details
        )
    elif NheadersFound < len(Headers):
        raise IOError(
            f"Not all of the Headers are in sampleTable, found only {NheadersFound}, something went wrong."
            + Error_details
        )

    files = sampleTable.loc[sample, Headers]

    if files.isnull().all():
        raise FileNotInSampleTableException(
            "The following files were not available for this sample in the SampleTable"
            + Error_details
        )

    elif files.isnull().any():
        raise IOError(
            f"Not all of the files are in sampleTable, something went wrong."
            + Error_details
        )

    return list(files)


def get_quality_controlled_reads(wildcards, include_se=False):
    """
    Gets quality controlled reads.
    R1 and R1 or se are returned as a dict.

    if the files are not in the sample tible impute default path produced with atlas.
    set

    """

    Fractions = MULTIFILE_FRACTIONS

    if config.get("interleaved_fastqs", False) and SKIP_QC:
        Fractions = ["se"]

    elif not include_se:
        # get only R1 and R2 or se
        Fractions = Fractions[: min(len(Fractions), 2)]

    try:
        QC_Headers = ["Reads_QC_" + f for f in Fractions]
        return get_files_from_sampleTable(wildcards.sample, QC_Headers)

    except FileNotInSampleTableException:
        # return files as named by atlas pipeline
        return expand(
            "QC/reads/{sample}_{fraction}.fastq.gz",
            fraction=Fractions,
            sample=wildcards.sample,
        )


def get_assembly(wildcards):
    """
    Returns Assembly file for a given sample.

    """

    Header = "Assembly"
    try:
        return get_files_from_sampleTable(wildcards.sample, Header)

    except FileNotInSampleTableException:
        # return files as named by atlas pipeline

        return "samples/{sample}/assembly/{sample}.fasta".format(sample=wildcards.sample)

logger.debug(f"\n{sampleTable}")
