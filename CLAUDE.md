# nested_RF_stimulus - Project Guide

## Experiment Overview

Bar flash stimuli presented at 8 orientations and 11 positions along the orientation axis to T4/T5 visual neurons in *Drosophila*. Recordings are intracellular voltage traces (whole-cell patch clamp, ~10 kHz sampling rate). Four experimental conditions are used: control vs TTL (pharmacological manipulation), each with stimulus ON and OFF epochs.

The 8 orientation columns correspond to bar orientations spaced around the circle. The "preferred-null" (PDND) axis is the orientation with the strongest spatial modulation. The "orthogonal direction" (OD) axis is 4 columns away (90 degrees offset).

## Repository Structure

```
src/
  analysis/
    helper/              - Utility functions (FWHM, circular variance, Gaussian fitting, colormaps)
    plotting/            - Visualization functions (see Plotting Conventions below)
    protocol2/           - Protocol 2 processing scripts
      pipeline/          - Data parsing functions (parse_bar_flash_data.m, etc.)
    results_analysis/    - Results analysis scripts
      analyze_bar_flash.m  - 1D RF bar flash analysis pipeline (PDND + OD)
    analyse_bar_DS/      - Bar direction selectivity analysis
  stimulus_generation/   - Stimulus pattern and position function generation
  protocol_generation/   - Experimental protocol scripts
  tests/                 - Test scripts
results/                 - Stimulus patterns, functions, and parameters
protocols/               - Experimental protocol definitions
```

## Git Branches

- `main` - Main branch
- `bar-flash-analysis` - 1D RF bar flash analysis pipeline (analyze_bar_flash.m + this CLAUDE.md)
- `setup-pharma` - Pharmacology experiment setup
- `process-summer-2025-data` - Summer 2025 data processing
- Other feature branches: `add-docstrings`, `test-AUC`, `visualisation_webpage_test`

## Data Structure

### Raw data location (not in repo)
`/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_flash_results/1DRF/`

Contains ~25 `.mat` files, one per recorded cell.

### .mat file contents
Each file contains (among other variables):
- `mean_slow`: 11 x 8 cell array (positions x orientations), each cell holds a time series (~15,800 samples at 10 kHz)
- `mean_fast`: Same structure for fast bar flash condition
- `data_slow`, `data_fast`: Individual repetition data (11 x 8 x 3 cell arrays, 3 reps)

### Important data characteristics
- **Voltage range**: Approximately -60 to -45 mV (all values are negative)
- **Time series lengths vary**: Most are ~15,801 samples, but some are 15,781 or 15,941
- **Some arrays are (1, N) shaped** instead of (N, 1) — the analysis script flattens these on load

### File naming convention
`bar_flash_results_YYYY_MM_DD_HH_MM_42F06_T4T5_{condition}.mat`
- Conditions determined by filename keywords: `control` vs `ttl`, and `_on` vs `_off`
- Groups into: `control_on`, `control_off`, `ttl_on`, `ttl_off`

## Key Analysis: analyze_bar_flash.m

Located at `src/analysis/results_analysis/analyze_bar_flash.m`.

### Pipeline
1. **Load & group** — Scans .mat files, groups by condition from filename keywords
2. **Column identification** — Finds preferred orientation per cell using two methods:
   - Method 1 (primary): Column with highest max voltage across all 11 positions
   - Method 2 (comparison): Column with largest max-min difference across position peaks
3. **Diagnostic plots** — Per-cell 2x4 subplot figures showing all 8 columns with peak markers and method annotations
4. **Reorder** — Circular shift so max-response position maps to position 6
5. **Peak alignment** — Centers voltage peaks at sample 10,000 in 20,000-length NaN arrays
6. **Averaging** — `nanmean` and SEM across cells per condition (for both PDND and OD axes)
7. **Per-condition figures** — 1x11 subplot figures with individual cell traces (gray) + mean (black)
8. **Comparison figures** — Control (black) vs TTL (red) mean + SEM shaded area, for both PDND and OD
9. **Summary report** — Text file with method comparison table, statistics, and OD column assignments

### Orthogonal direction (OD) analysis
Once the PDND column is found, the OD column is computed as `mod(best_col - 1 + 4, 8) + 1`. The full reorder → align → average pipeline runs in parallel for both axes.

### Output files
- `bar_flash_analysis_PDND_{condition}.fig/.png` — Per-condition PDND figures
- `column_diagnostic_{condition}_cell{NN}.fig/.png` — Column selection diagnostic per cell
- `bar_flash_comparison_PDND_on.fig/.png` — Control vs TTL comparison, ON stimulus, PDND axis
- `bar_flash_comparison_PDND_off.fig/.png` — Control vs TTL comparison, OFF stimulus, PDND axis
- `bar_flash_comparison_OD_on.fig/.png` — Control vs TTL comparison, ON stimulus, OD axis
- `bar_flash_comparison_OD_off.fig/.png` — Control vs TTL comparison, OFF stimulus, OD axis
- `analysis_summary.txt` — Summary report
- `processed_data.mat` — All aligned and averaged data (PDND + OD)

### Key parameters
```matlab
num_positions = 11;      % Rows in mean_slow
num_orientations = 8;    % Columns in mean_slow
center_position = 6;     % Target position for max response after circular shift
aligned_length = 20000;  % Length of peak-aligned NaN arrays
peak_center = 10000;     % Sample index where peak is centered
od_offset = 4;           % Column offset for orthogonal direction
```

## Known Pitfalls and Bugs Found

### Negative voltage initialization
When finding max voltage across time series, arrays **must** be initialized to `-Inf`, not `zeros`. Since all voltage values are negative (~-60 to -45 mV), initializing to 0 means `max(0, max(ts))` always returns 0, causing the first column to be selected by default regardless of actual data. This was a bug in the initial implementation of `find_preferred_column_method1` and `reorder_to_center_max`.

### Method 2 is not affected by the same bug
`find_preferred_column_method2` computes differences (`max - min`), which are positive values, so initializing `pos_maxes` to zeros and then overwriting with `max(ts)` works correctly since the values get fully replaced.

## Plotting Conventions

Based on existing scripts in `src/analysis/plotting/`:

- **Layout**: `tiledlayout`/`nexttile` for subplot grids
- **Axes style**: `box off`, `ax.TickDir = 'out'`, `ax.TickLength = [0.04 0.04]`
- **Y-axis**: Visible only on first subplot; `ylabel('Voltage (mV)')`. Others: `ax.YAxis.Visible = 'off'`
- **X-axis**: Often hidden or labeled with time in seconds; `xlabel('Time (s)')`
- **Colors**: Control/baseline in black (`'k'`), experimental/TTL in red (`'r'`)
- **Individual traces**: Gray `[0.8 0.8 0.8]`, thin `LineWidth` 0.5–0.7
- **Mean traces**: `LineWidth` 1.5–2
- **SEM**: Shaded `fill()` regions with `FaceAlpha` 0.4, gray for control, light red `[1 0.8 0.8]` for TTL
- **Legends**: `'Box', 'off'`, small font, placed in first subplot only
- **Export**: Both `.fig` (MATLAB) and `.png` (300 dpi via `exportgraphics`); some scripts also export PDF with `'ContentType', 'vector'`
- **Figure positions**: Wide horizontal for 1xN subplot layouts (e.g., `[50 200 2000 400]`)

## Coding Conventions

- MATLAB scripts and functions (not Python)
- Cell arrays for variable-length time series
- `nanmean`/`nanstd` for handling NaN-padded arrays
- Local functions at bottom of scripts (not separate files)
- `circshift` for circular reordering of position arrays
- Condition grouping via filename string matching (`contains()`)
- `fprintf` for console progress reporting during long pipelines
- `-v7.3` flag when saving large `.mat` files

## Data Processing Pipeline (upstream)

The raw electrophysiology data is processed by `parse_bar_flash_data.m` (in `src/analysis/protocol2/pipeline/`):
1. Frame data sampled at 10 kHz identifies stimulus timing
2. 3 repetitions per stimulus condition
3. Bar flash epochs extracted using gap detection (3s grey screen intervals)
4. Per-position, per-orientation time series extracted and averaged across reps
5. Output: `mean_slow` and `mean_fast` cell arrays stored in per-cell `.mat` files
