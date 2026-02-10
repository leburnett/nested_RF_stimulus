# nested_RF_stimulus - Project Guide

## Experiment Overview

Bar flash stimuli presented at 8 orientations and 11 positions along the orientation axis to T4/T5 visual neurons in *Drosophila*. Recordings are intracellular voltage traces (whole-cell patch clamp). Four experimental conditions are used: control vs TTL (pharmacological manipulation), each with stimulus ON and OFF epochs.

## Repository Structure

```
src/
  analysis/
    helper/              - Utility functions (FWHM, circular variance, etc.)
    plotting/            - Visualization functions
    protocol2/           - Protocol 2 processing scripts
      pipeline/          - Data parsing functions (parse_bar_flash_data.m, etc.)
    results_analysis/    - Results analysis scripts
      analyze_bar_flash.m  - 1D RF bar flash analysis pipeline
  stimulus_generation/   - Stimulus pattern and position function generation
  protocol_generation/   - Experimental protocol scripts
results/                 - Stimulus patterns, functions, and parameters
protocols/               - Experimental protocol definitions
```

## Data Structure

### Raw data location (not in repo)
`/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_flash_results/1DRF/`

### .mat file contents
Each file contains (among other variables):
- `mean_slow`: 11 x 8 cell array (positions x orientations), each cell holds a time series (~15,800 samples at 10 kHz)
- `mean_fast`: Same structure for fast bar flash condition
- `data_slow`, `data_fast`: Individual repetition data (11 x 8 x 3 cell arrays)

### File naming convention
`bar_flash_results_YYYY_MM_DD_HH_MM_42F06_T4T5_{condition}.mat`
- Conditions: `control_on`, `control_off`, `ttl_on`, `ttl_off`

## Key Analysis: analyze_bar_flash.m

Located at `src/analysis/results_analysis/analyze_bar_flash.m`.

### Pipeline
1. **Load & group** - Scans .mat files, groups by condition from filename
2. **Column identification** - Finds preferred orientation per cell:
   - Method 1 (primary): Column with highest max voltage across all positions
   - Method 2 (comparison): Column with largest max-min difference across positions
3. **Reorder** - Circular shift so max-response position maps to position 6
4. **Peak alignment** - Centers voltage peaks at sample 10,000 in 20,000-length NaN arrays
5. **Averaging** - nanmean and SEM across cells per condition
6. **Visualization** - 1x11 subplot figures (individual traces + mean overlay)
7. **Summary report** - Text file with method comparison and statistics

### Outputs
- `bar_flash_analysis_{condition}.fig` / `.png` - Figures per condition
- `analysis_summary.txt` - Summary report
- `processed_data.mat` - All aligned and averaged data

## Coding Conventions
- MATLAB scripts and functions
- Cell arrays for variable-length time series
- `tiledlayout`/`nexttile` for subplot grids
- `nanmean`/`nanstd` for handling NaN-padded arrays
- Local functions at bottom of scripts
