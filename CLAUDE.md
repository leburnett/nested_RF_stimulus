# CLAUDE.md — nested_RF_stimulus

Documentation of the stimulus generation, protocol structure, and analysis mapping chain for Protocol 2 (1DRF experiments). Written to prevent re-discovery of the relationships described below.

---

## Project Overview

MATLAB neuroscience analysis pipeline for Drosophila visual neuron whole-cell patch-clamp recordings. Stimuli are presented on a G4 LED arena. Data is sampled at 10 kHz (voltage + frame channel) via TDMS files.

**Protocol 2** consists of:
1. Small-field flashes (RF mapping)
2. Bar sweeps: 8 bar orientations x 2 directions (forward + flip) x 3 speeds (28/56/168 dps) = 48 conditions
3. Bar flashes: 8 orientations x 11 positions x 3 repetitions

---

## Full-Field Bar Source Patterns

**Location:** `results/patterns/protocol2/full_field_bars4/`

16 files (sorted alphabetically by filename, which is the order MATLAB's `dir('*.mat')` returns):

| Index | Filename | Polarity | Rotation Angle |
|-------|----------|----------|----------------|
| 0001 | `0001_4pix_bar_15_4_pi_39shift...` | ON (px=15) | pi |
| 0002 | `0002_4pix_bar_15_4_1-8pi_44shift...` | ON | pi/8 |
| 0003 | `0003_4pix_bar_15_4_2-8pi_44shift...` | ON | 2pi/8 |
| 0004 | `0004_4pix_bar_15_4_3-8pi_44shift...` | ON | 3pi/8 |
| 0005 | `0005_4pix_bar_15_4_4-8pi_44shift...` | ON | 4pi/8 |
| 0006 | `0006_4pix_bar_15_4_5-8pi_44shift...` | ON | 5pi/8 |
| 0007 | `0007_4pix_bar_15_4_6-8pi_44shift...` | ON | 6pi/8 |
| 0008 | `0008_4pix_bar_15_4_7-8pi_44shift...` | ON | 7pi/8 |
| 0009 | `0009_4pix_bar_0_4_pi_44shift...` | OFF (px=0) | pi |
| 0010 | `0010_4pix_bar_0_4_1-8pi_44shift...` | OFF | pi/8 |
| 0011 | `0011_4pix_bar_0_4_2-8pi_44shift...` | OFF | 2pi/8 |
| 0012 | `0012_4pix_bar_0_4_3-8pi_44shift...` | OFF | 3pi/8 |
| 0013 | `0013_4pix_bar_0_4_4-8pi_44shift...` | OFF | 4pi/8 |
| 0014 | `0014_4pix_bar_0_4_5-8pi_44shift...` | OFF | 5pi/8 |
| 0015 | `0015_4pix_bar_15_4_6-8pi_44shift...` | ON* | 6pi/8 |
| 0016 | `0016_4pix_bar_15_4_7-8pi_44shift...` | ON* | 7pi/8 |

*Note: 0015-0016 have px=15 (ON) in the filename despite being in the OFF index range (9-16). This may be a naming inconsistency.

**Critical:** The rotation angles (pi, pi/8, 2pi/8, ...) are rotation angles applied to a **vertical bar** by a custom lab rotation script. They are NOT `imrotate` angles and the rotation convention (CW vs CCW) is unknown. The actual bar orientation in arena coordinates cannot be deduced from the filename alone — it must be verified visually.

Each pattern is a `48 x 192 x 288` array (height x width x frames), gs_val=4 (pixel values 0-15).

---

## Experiment Pattern Generation

Orchestrated by `src/protocol_generation/generate_protocol2.m`:

### Step 1: Flash patterns (pattern files 0001-0002)
- Created by `generate_flash_stimulus_xy.m`
- Small-field ON/OFF flashes for RF mapping

### Step 2: Bar sweep patterns (pattern files 0003-0010)
- Created by `src/stimulus_generation/pattern/generate_bar_stimulus_xy.m`
- Loads full-field bars via `bar_patts = dir('*.mat')` (alphabetical order)
- For ON cells: `patt_ids = 1:8` → loads full-field indices 0001-0008
- For each, calls `generate_bar_patterns_xy.m` which crops the full-field pattern to a `px_crop` square centred on the cell's [x, y] RF coordinate
- **Mapping: experiment pattern 0003 = full-field 0001 (pi), 0004 = full-field 0002 (pi/8), ..., 0010 = full-field 0008 (7pi/8)**

### Step 3: Bar flash pattern (pattern file 0011)
- Created by `src/stimulus_generation/pattern/generate_bar_flash_stimulus_xy.m`
- Same `patt_ids = 1:8`, same `dir()` alphabetical order as sweep patterns
- Iterates p=1 through p=8, filling frames sequentially:

| Loop p | Full-field loaded | Rotation | Output frames (0-indexed) | Pats indices (+1 for grey bg) |
|--------|-------------------|----------|---------------------------|-------------------------------|
| 1 | 0001 | pi | 1-11 | 2-12 |
| 2 | 0002 | pi/8 | 12-22 | 13-23 |
| 3 | 0003 | 2pi/8 | 23-33 | 24-34 |
| 4 | 0004 | 3pi/8 | 34-44 | 35-45 |
| 5 | 0005 | 4pi/8 | 45-55 | 46-56 |
| 6 | 0006 | 5pi/8 | 56-66 | 57-67 |
| 7 | 0007 | 6pi/8 | 67-77 | 68-78 |
| 8 | 0008 | 7pi/8 | 78-88 | 79-89 |

- `centre_frames = [30, 33, 33, 34, 34, 34, 34, 35, ...]` — the frame in each full-field pattern where the bar is centred on the arena
- Frame extraction: `frames_from_patt = cf - (n_flank*2) : 2 : cf + (n_flank*2)` — takes every other frame to get 11 positions spanning the bar's sweep range
- Frame 1 of the output pattern (Pats index 1) is a grey background

**Key result: bar flash column 1 = pi rotation, column 2 = pi/8 rotation, ..., column 8 = 7pi/8 rotation. This is the same sequence as sweep patterns 0003-0010.**

---

## Protocol Order and Data Row Mapping

### Protocol order (from `create_protocol2.m`, experiment-era version commit 0181053)

```
flash_patt = [1, 2]    % pattern files 0001, 0002
bar_patt = 3:10        % pattern files 0003-0010 (8 orientations)
n_speeds = 3           % 28, 56, 168 dps
```

**Warning:** The current branch version of `create_protocol2.m` has `flash_patt = 1`, `bar_patt = 2:10`, `n_speeds = 5`. This does NOT match the 1DRF data which was collected with the 3-speed version.

### Protocol sequence for bar sweeps

For each speed, patterns 3-10 are played with alternating position functions:
- Function 3 = forward (frames 11→62)
- Function 4 = flip (frames 62→11)

```
Entry 5:  pattern 3, func 3  (forward)
Entry 6:  pattern 3, func 4  (flip)
Entry 7:  pattern 4, func 3  (forward)
Entry 8:  pattern 4, func 4  (flip)
...
Entry 19: pattern 10, func 3 (forward)
Entry 20: pattern 10, func 4 (flip)
```

### Data row assignment (from `parse_bar_data.m`)

`parse_bar_data` extracts voltage timeseries sequentially from the recording. Within each speed block, data rows 1-16 correspond to:

| Data Row | Pattern File | Direction | Rotation |
|----------|-------------|-----------|----------|
| 1 | 0003 | Forward | pi |
| 2 | 0003 | Flip | pi |
| 3 | 0004 | Forward | pi/8 |
| 4 | 0004 | Flip | pi/8 |
| 5 | 0005 | Forward | 2pi/8 |
| 6 | 0005 | Flip | 2pi/8 |
| ... | ... | ... | ... |
| 15 | 0010 | Forward | 7pi/8 |
| 16 | 0010 | Flip | 7pi/8 |

**Rule: odd data rows = forward, even data rows = flip. Pattern file index = ceil(data_row / 2) + 2.**

---

## Polar Plot Mapping

### `plot_order` (in `plot_timeseries_polar_bars.m`)

```matlab
plot_order = [1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 10, 12, 14, 16];
angls = linspace(0, 2*pi, 17); angls = angls(1:16);
```

Maps **angular index** (subplot position on polar plot) → **data row**:

| Ang Idx | Angle (deg) | Data Row | Direction | Pattern |
|---------|-------------|----------|-----------|---------|
| 1 | 0.0 | 1 | Forward | 0003 (pi) |
| 2 | 22.5 | 3 | Forward | 0004 (pi/8) |
| 3 | 45.0 | 5 | Forward | 0005 (2pi/8) |
| 4 | 67.5 | 7 | Forward | 0006 (3pi/8) |
| 5 | 90.0 | 9 | Forward | 0007 (4pi/8) |
| 6 | 112.5 | 11 | Forward | 0008 (5pi/8) |
| 7 | 135.0 | 13 | Forward | 0009 (6pi/8) |
| 8 | 157.5 | 15 | Forward | 0010 (7pi/8) |
| 9 | 180.0 | 2 | Flip | 0003 (pi) |
| 10 | 202.5 | 4 | Flip | 0004 (pi/8) |
| 11 | 225.0 | 6 | Flip | 0005 (2pi/8) |
| 12 | 247.5 | 8 | Flip | 0006 (3pi/8) |
| 13 | 270.0 | 10 | Flip | 0007 (4pi/8) |
| 14 | 292.5 | 12 | Flip | 0008 (5pi/8) |
| 15 | 315.0 | 14 | Flip | 0009 (6pi/8) |
| 16 | 337.5 | 16 | Flip | 0010 (7pi/8) |

Angular indices 1-8 (0°-157.5°) = forward directions. Angular indices 9-16 (180°-337.5°) = flip (null) directions.

**Note from source code comment:** *"Dark bars - two directions are mixed up. Occurred when making patterns with bkg4..."* — this may indicate `plot_order` was adjusted to account for a known issue but the adjustment may not be fully correct.

### Response magnitude

```matlab
max_v(i, sp) = abs(diff([prctile(d, 98), mean_before]));
```

Measures the absolute magnitude of the 98th percentile depolarisation relative to the pre-stimulus baseline. PD is the direction/speed combination with the highest `max_v` value.

---

## PD → Pattern → Bar Flash Mapping Chain

```
1. PD detection:      [~, idx] = max(max_v(:))
                      [pd_dir_idx, pd_speed_idx] = ind2sub(size(max_v), idx)
                      pd_angle_rad = angls(pd_dir_idx)

2. Data row:          pd_data_row = plot_order(pd_dir_idx)

3. Pattern file:      pattern_file_idx = ceil(pd_data_row / 2) + 2
                      is_flip = mod(pd_data_row, 2) == 0

4. Bar flash column:  bf_col = ceil(pd_data_row / 2)
```

**Why the direct mapping works:** Sweep patterns 0003-0010 and bar flash columns 1-8 are both generated by iterating `patt_ids = 1:8` over the same alphabetically-ordered full-field bar files. So sweep pattern `i+2` and bar flash column `i` both use full-field bar `i`, which has the same rotation angle. No trigonometry is needed.

**Previous bug (now fixed):** The old `map_PD_to_bar_flash_column` treated the rotation angles as absolute bar orientations and used `bar_orient = mod(pd_angle_rad - pi/2, pi)` to find the closest match. This produced **orthogonal** bar orientations because the rotation angles are not absolute orientations.

---

## Bar Flash Data Structure

`parse_bar_flash_data.m` returns `data_slow{position, orientation_col, rep}` — an 11x8x3 cell array.

- **Dimension 1:** 11 positions (bar locations along the sweep axis)
- **Dimension 2:** 8 orientation columns (matching the 8 rotation angles)
- **Dimension 3:** 3 repetitions

Frame number formula: `frame = (col - 1) * 11 + pos + 1`

The `+1` offset accounts for frame 1 being a grey background in pattern 0011.

`parse_bar_flash_data` uses `flash_frame_num = max(f_data(st:nd))` to determine which cell to store the voltage data in. The position function randomises the frame order during the experiment.

---

## Known Pitfalls

1. **Voltage values are all negative** (~-60 to -45 mV). Initialise max-finding to `-Inf`, not zeros.
2. **`load_protocol2_data` uses `cd` internally.** Always use `onCleanup` to restore the working directory.
3. **`plot_timeseries_polar_bars` expects `params.date`, `params.time`, `params.strain`** but `load_protocol2_data` does not set these. They must be added manually:
   ```matlab
   params.date = date_str;
   params.time = time_str;
   loaded_exp = load(fullfile(exp_folder, 'currentExp.mat'));
   params.strain = loaded_exp.metadata.Strain;
   ```
4. **MATLAB `imagesc` y-axis inversion.** Default is `'reverse'` (row 1 at top). The physical arena may use the opposite convention. Stimulus display figures may show vertical motion inverted. Use `set(gca, 'YDir', 'normal')` to test the alternative.
5. **`create_protocol2.m` on the current branch does not match 1DRF data.** The experiment-era version (commit `0181053`) had `flash_patt = [1, 2]`, `bar_patt = 3:10`, `n_speeds = 3`.
6. **The `rotation_angles` array in `map_PD_to_bar_flash_column.m` is for display labels only.** The actual mapping uses `bf_col = ceil(pd_data_row / 2)` which is independent of angle values.

---

## Active Investigation: Direction Alignment

The sweep stimulus in debug figures appears to show the **null direction** (not PD) for at least some cells. Possible causes:

1. **Y-axis display issue:** MATLAB's `imagesc` default `'reverse'` y-direction may invert vertical motion relative to the arena. The `debug_PD_mapping.m` script accepts a `y_dir` parameter (`'reverse'` or `'normal'`) to toggle this.

2. **`plot_order` may have incorrect angular assignments:** The source comment about "dark bars - two directions are mixed up" suggests a known issue. If the forward/flip angular assignments are swapped, the detected PD would be 180° off.

3. **Interactive debugging:** `debug_PD_mapping.m` Step 13 provides a 4x4 grid of all 16 directions and an interactive prompt to manually select which direction you believe is the PD. The script then traces the full mapping chain from that selection.

---

## Key Files

### Stimulus Generation
| File | Description |
|------|-------------|
| `src/protocol_generation/generate_protocol2.m` | Orchestrates full protocol 2 setup |
| `src/protocol_generation/create_protocol2.m` | Creates `currentExp.mat` with pattern/func orders |
| `src/stimulus_generation/pattern/generate_bar_stimulus_xy.m` | Creates cropped bar sweep patterns (0003-0010) |
| `src/stimulus_generation/pattern/generate_bar_patterns_xy.m` | Crops single full-field bar to experiment coordinates |
| `src/stimulus_generation/pattern/generate_bar_flash_stimulus_xy.m` | Creates bar flash pattern (0011) from full-field bars |
| `src/stimulus_generation/posn_fn/generate_bar_pos_fns.m` | Creates forward/flip position functions per speed |
| `src/stimulus_generation/posn_fn/generate_bar_flash_pos_fns.m` | Creates randomised bar flash position functions |

### Analysis Pipeline
| File | Description |
|------|-------------|
| `src/analysis/protocol2/pipeline/load_protocol2_data.m` | Loads TDMS data and params (uses `cd` internally) |
| `src/analysis/protocol2/pipeline/parse_bar_data.m` | Parses bar sweep voltage data → 48x4 cell |
| `src/analysis/protocol2/pipeline/parse_bar_flash_data.m` | Parses bar flash data → 11x8x3 cell |
| `src/analysis/plotting/plot_timeseries_polar_bars.m` | Polar plot + compute max_v/min_v |
| `src/analysis/helper/map_PD_to_bar_flash_column.m` | Maps PD data row → bar flash column (direct mapping) |
| `src/analysis/plotting/plot_bar_flash_PD_positions.m` | Bar flash spatial profile plots (raw + baseline-subtracted) |
| `src/analysis/plotting/create_bar_sweep_gif.m` | Animated GIF of bar sweep stimulus |
| `src/analysis/protocol2/process_1DRF_single.m` | Single experiment processing pipeline |
| `src/analysis/protocol2/process_1DRF.m` | Batch processing for all 1DRF experiments |
| `src/analysis/debug_PD_mapping.m` | Interactive diagnostic for the full mapping chain |

### Data Paths
| Path | Description |
|------|-------------|
| `results/patterns/protocol2/full_field_bars4/` | Source full-field bar patterns (16 files) |
| `/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/1DRF/` | Experiment data (26 cells) |
| `/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/1DRF/` | Analysis results |

### Branch Structure
| Branch | Description |
|--------|-------------|
| `main` | Standard protocol |
| `setup-pharma` | Pharmacology variant (5 speeds) |
| `bar-flash-analysis` | Bar flash analysis development |
| `feature/process-1DRF` | Batch 1DRF processing pipeline (current work) |
