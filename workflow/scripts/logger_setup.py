import sys
import logging
from snakemake.utils import logger
from naive_atlas.color_logger import ColorFormatter, logging_format, datefmt

def configure_snakemake_logger(config):
    """Configures Snakemake's internal loggers with custom color formatting."""
    
    # Create clean console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(ColorFormatter(logging_format, datefmt=datefmt))
    
    # Target every logging namespace Snakemake uses
    target_loggers = [
        logging.getLogger(),                   # The system root logger
        logging.getLogger("snakemake"),        # Base snakemake framework logger
        logging.getLogger("snakemake.logging") # Core engine logging module
    ]
    
    # If the imported wrapper object contains a nested logger, include it
    if hasattr(logger, "logger") and logger.logger is not None:
        target_loggers.append(logger.logger)
    
    for target in target_loggers:
        # Remove all default handlers so they can't print raw text
        for handler in list(target.handlers):
            # target.removeHandler(handler)
            # target.addHandler(console_handler)
            target.propagate = False # Prevent messages bubbling up to double-print
    
    # Force ignore DEBUG chatter even when global debug is active
    logging.getLogger("asyncio").setLevel(logging.WARNING)
    logging.getLogger("smart_open.smart_open_lib").setLevel(logging.WARNING)
    logging.getLogger("py.warnings").setLevel(logging.WARNING)
    
    # Set appropriate logging thresholds based on your debug flag
    is_debug = str(config.get("debug", "False")).lower() in ("true", "1", "yes")
    
    if is_debug:
        for target in target_loggers:
            target.setLevel(logging.DEBUG)
        console_handler.setLevel(logging.DEBUG)
        logger.stdout = True
    else:
        for target in target_loggers:
            target.setLevel(logging.INFO)
        console_handler.setLevel(logging.INFO)
        logger.stdout = False

