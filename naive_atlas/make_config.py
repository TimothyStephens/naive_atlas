from .default_values import *
from snakemake.utils import update_config as snakemake_update_config
from snakemake.common.configfile import load_configfile
import tempfile
import sys
import os
import multiprocessing
import logging

logger = logging.getLogger(__file__)



def make_config(
    database_dir,
    config="config.yaml",
):
    """
    Reads template config file with comments from ../workflow/config/template_config.yaml
    updates it by the parameters provided.

    Args:
        config (str): output file path for yaml
        database_dir (str): location of downloaded databases
    """

    from ruamel.yaml import YAML  # used for yaml reading with comments

    yaml = YAML()

    template_conf_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "workflow/../config/template_config.yaml",
    )

    with open(template_conf_file) as template_config:
        conf = yaml.load(template_config)

    if os.path.exists(config):
        logger.warning(
            f"Config file {config} already exists, I didn't dare to overwrite it. continue..."
        )
    else:
        with open(config, "w") as f:
            yaml.dump(conf, f)
        logger.info(
            "Configuration file written to %s\n"
            "        You may want to edit it using any text editor." % config
        )


