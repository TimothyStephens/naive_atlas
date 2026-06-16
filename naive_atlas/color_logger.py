import sys
import os
import logging
import traceback

# ==============================================================================
# 1. Dynamically Safely Fetch Snakemake's Native Colors
# ==============================================================================
prefix = ""

# Fetch Snakemake's live engine logger
snakemake_base_logger = logging.getLogger("snakemake")

# Safe lookups: Try to grab Snakemake's live color codes from its active handlers.
# If they aren't available, automatically fall back to standard ANSI values.
def get_native_color(level, default_ansi):
    if hasattr(snakemake_base_logger, "logger") and snakemake_base_logger.logger:
        for handler in snakemake_base_logger.logger.handlers:
            if hasattr(handler, "colors") and level in handler.colors:
                return handler.colors[level]
    return default_ansi

# Bind to Snakemake's palette or use classic clean ANSI fallbacks
color_debug  = get_native_color(logging.DEBUG, "\033[0;35m")   # Purple / Dark Gray
color_info   = get_native_color(logging.INFO, "\033[0;32m")    # Green
color_warn   = get_native_color(logging.WARNING, "\033[0;33m") # Yellow
color_error  = get_native_color(logging.ERROR, "\033[0;31m")   # Red
color_reset  = "\033[0m"

logging_format = "[%(asctime)s - %(name)s - %(levelname)s]: %(message)s"
datefmt = "%Y-%m-%d %H:%M"

# ==============================================================================
# 2. Custom Color Formatter
# ==============================================================================
class ColorFormatter(logging.Formatter):
    def __init__(self, format=logging_format, datefmt=datefmt):
        super().__init__(fmt=format, datefmt=datefmt)
        
        self.FORMATS = {
            logging.DEBUG: prefix + color_debug + format + color_reset,
            logging.INFO: prefix + color_info + format + color_reset,
            logging.WARNING: prefix + color_warn + format + color_reset,
            logging.ERROR: prefix + color_error + format + color_reset,
            logging.CRITICAL: prefix + color_error + format + color_reset,
        }

    def format(self, record):
        if record.name.startswith("snakemake"):
            record.name = "snakemake"
        elif record.name.startswith("naive_atlas"):
            record.name = "naive_atlas"

        log_fmt = self.FORMATS.get(record.levelno, prefix + self._fmt + color_reset)
        formatter = logging.Formatter(log_fmt, datefmt=self.datefmt)
        return formatter.format(record)

# ==============================================================================
# 3. Base Logger Configuration
# ==============================================================================
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setLevel(logging.DEBUG)  
consoleHandler.setFormatter(ColorFormatter(logging_format, datefmt=datefmt))

logging.basicConfig(
    level=logging.INFO,  
    handlers=[consoleHandler],
)

logger = logging.getLogger("naive_atlas")
logger.addHandler(consoleHandler)
logger.propagate = False  

logging.captureWarnings(True)

# ==============================================================================
# 4. Uncaught Exception Hook Handler
# ==============================================================================
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    error_msg = "".join(
        traceback.format_exception(exc_type, exc_value, exc_traceback)
    )
    logger.critical(f"Uncaught exception detected:\n{error_msg}")

sys.excepthook = handle_exception

