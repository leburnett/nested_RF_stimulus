# Nested Receptive Field Protocol Visualization

This folder contains a web-based visualization explaining the two-step nested receptive field protocol.

## Quick Start

1. **Generate the GIFs** (requires MATLAB):
   ```matlab
   cd('/Users/burnettl/Documents/GitHub/nested_RF_stimulus/src/visualization')
   generate_all_gifs
   ```

2. **Add the required images** to the `assets/` folder:
   - `arena_setup.png` - The G4 arena experimental setup image
   - `analysis_dark.png` - Response grid for dark (OFF) flashes
   - `analysis_light.png` - Response grid for light (ON) flashes

3. **Open the visualization**:
   - Open `index.html` in Chrome (or any modern browser)
   - File > Open File > select `index.html`

## File Structure

```
visualization/
├── index.html              # Main webpage
├── styles.css              # CSS styling
├── generate_stimulus_gif.m # Updated GIF generator function
├── generate_all_gifs.m     # Script to generate all required GIFs
├── README.md               # This file
└── assets/                 # Images and GIFs
    ├── arena_setup.png     # (you provide)
    ├── analysis_dark.png   # (you provide)
    ├── analysis_light.png  # (you provide)
    ├── protocol1_12px.gif  # (generated)
    ├── protocol1_6px.gif   # (generated)
    ├── protocol2_region.gif# (generated)
    ├── protocol2_4px.gif   # (generated)
    ├── protocol2_6px.gif   # (generated)
    └── protocol2_bars.gif  # (generated)
```

## Customization

### Changing the demo peak location

Edit `generate_all_gifs.m` and modify these variables:
```matlab
demo_peak_x = 64;  % X coordinate of demo peak
demo_peak_y = 24;  % Y coordinate of demo peak
```

### Adjusting GIF speed/length

Each GIF has options you can adjust in `generate_all_gifs.m`:
- `fps_initial` - Starting playback speed
- `fps_fast` - Accelerated playback speed
- `switch_time_step` - When to switch from slow to fast
- `max_time_steps` - How much of the protocol to show
- `time_step_skip` - Subsampling factor (higher = smaller file)

## Notes

- The Protocol 2 flash GIFs currently use Protocol 1 patterns as proxies (cropped to the relevant region). For accurate Protocol 2 visualizations, first generate actual Protocol 2 patterns using `generate_protocol2.m`.

- The visualization uses a dark theme with green coloring to match the appearance of the real G4 arena LEDs.
