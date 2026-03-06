"""Constants for the nested RF stimulus ephys dashboard."""

# Default directory for preprocessed dashboard images
DEFAULT_IMAGES_DIR = (
    "/Users/burnettl/Documents/Projects/nested_RF_stimulus/dashboard_images"
)

# Experiment log path
EXPERIMENT_LOG_PATH = (
    "/Users/burnettl/Documents/GitHub/nested_RF_stimulus/"
    "reiser-lab-ephys-experiment-log.xlsx"
)

# Dashboard port (8051 to avoid conflict with freely-walking on 8050)
PORT = 8051

# Project type filters for scoping
PROJECT_TYPES = ["Summer2025", "Autumn2025"]

# Strain to exclude
EXCLUDED_STRAINS = ["42F06_T4T5_test"]

# Image types and their display names
IMAGE_TYPES = {
    "gridplot_1": "Grid Plot 1",
    "gridplot_2": "Grid Plot 2",
    "gridplot_3": "Grid Plot 3",
    "gridplot_4": "Grid Plot 4",
    "bar_polar": "Bar Sweep Polar",
    "flash_rf": "Square Flash RF Timeseries",
    "flash_rf_heatmap": "Flash RF Heatmap",
    "bar_flash": "Bar Flash Timeseries",
}

# Number of gridplot sub-figures per cell
N_GRIDPLOTS = 4

# Cache max-age for static images (seconds) — 24 hours
IMAGE_CACHE_MAX_AGE = 86400
