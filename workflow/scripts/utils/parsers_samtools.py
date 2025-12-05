import pandas as pd


def parse_samtools_stats_file(stats_file):
    """
    parses a samtools stats file (paired or single end)
    returns a dixt with parsed values
    """

    parsed = {}

    with open(stats_file) as f:
        for line in f:
            if line.startswith('SN'):
                try:
                    lsplit = line.rstrip().split("\t")
                    group, key, value = lsplit[:3]

                    try:
                        value = int(value)
                    except ValueError:
                        value = float(value)

                    parsed[key.strip(":")] = value

                except Exception as e:
                    raise Exception(
                        f"Error parsing line:\n{line}\n in log file {stats_file}"
                    ) from e

        return parsed


