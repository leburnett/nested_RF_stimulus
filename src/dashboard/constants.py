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
    "polar_arrow": "Polar with Arrow",
    "flash_rf_4px": "Flash RF Timeseries (4px)",
    "flash_heatmap_4px": "Flash RF Heatmap (4px)",
    "flash_rf_6px": "Flash RF Timeseries (6px)",
    "flash_heatmap_6px": "Flash RF Heatmap (6px)",
    "gaussian_rf": "Gaussian RF Contours",
    "bar_flash_slow": "Bar Flash (80ms)",
    "bar_flash_fast": "Bar Flash (14ms)",
}

# Number of gridplot sub-figures per cell
N_GRIDPLOTS = 4

# Cache max-age for static images (seconds) — 24 hours
IMAGE_CACHE_MAX_AGE = 86400
