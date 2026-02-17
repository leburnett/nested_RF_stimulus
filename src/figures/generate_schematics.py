#!/usr/bin/env python3
"""Generate schematic figures of LED arena stimulus protocols.

Generates 8 figures showing Protocol 1 and Protocol 2 stimulus configurations
for the nested RF mapping system. Each figure is saved as both SVG and PDF.

Usage:
    python src/figures/generate_schematics.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

# ============================================================================
# Constants
# ============================================================================

# Arena dimensions
ARENA_W = 192
ARENA_H = 48

# Valid display area (1-indexed, inclusive)
SCREEN_COL_START = 17
SCREEN_COL_END = 192
SCREEN_ROW_START = 1
SCREEN_ROW_END = 48

# Protocol 1 LHS display region
P1_COL_START = 17
P1_COL_END = 112
P1_ROW_START = 1
P1_ROW_END = 48

# Pixel intensity values (gs_val=4, range 0-15)
BKG_INTENSITY = 4
OFF_INTENSITY = 0
ON_INTENSITY = 15

# Protocol 2 parameters
P2_CROP_SIZE = 30
P2_FLASH_4PX = 4
P2_FLASH_6PX = 6
P2_OVERLAP = 0.5
BAR_WIDTH = 4
N_FLANK = 5
N_BAR_POSITIONS = 11
BAR_POSITION_STEP = 2

# Example RF center for Protocol 2 schematics (center of LHS region)
RF_CENTER_X = 65
RF_CENTER_Y = 24

# Grid line style
GRID_COLOR = 'lightgrey'
GRID_LW = 0.5
OUTLINE_COLOR = 'lightgrey'
OUTLINE_LW = 1.0

# Output
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
OUTPUT_DIR = os.path.join(REPO_ROOT, 'results', 'images')


# ============================================================================
# Helper functions
# ============================================================================

def intensity_to_rgb(val):
    """Map pixel intensity (0-15) to RGB. Black to bright green."""
    g = val / 15.0
    return (0.0, g, 0.0)


def make_arena_cmap():
    """Create a 16-entry colormap mapping intensity 0-15 to black-green."""
    colors = [intensity_to_rgb(i) for i in range(16)]
    return ListedColormap(colors)


ARENA_CMAP = make_arena_cmap()


def create_arena_array(fill_value=0):
    """Return a 48x192 array filled with a given intensity value."""
    return np.full((ARENA_H, ARENA_W), float(fill_value))


def fill_region(arr, row_start, row_end, col_start, col_end, value):
    """Fill a rectangular region using 1-indexed inclusive coordinates."""
    arr[row_start - 1:row_end, col_start - 1:col_end] = value


def centered_square(x, y, square_size):
    """Port of centeredSquare.m. Returns 1-indexed inclusive bounds."""
    half = square_size / 2
    row_start = int(max(SCREEN_ROW_START,
                        min(y - half + 1, SCREEN_ROW_END - square_size + 1)))
    row_end = row_start + square_size - 1
    col_start = int(max(SCREEN_COL_START,
                        min(x - half + 1, SCREEN_COL_END - square_size + 1)))
    col_end = col_start + square_size - 1
    return row_start, row_end, col_start, col_end


def compute_flash_positions(region_start, region_end, flash_size, overlap):
    """Compute flash edge start positions (1-indexed).

    Ports the grid calculation from generate_flash_pattern.m lines 78-85.
    """
    step = int(flash_size * (1 - overlap))
    return list(range(region_start, region_end - flash_size + 2, step))


def create_figure():
    """Create a figure sized for the 192x48 arena (4:1 aspect)."""
    fig_width = 10
    fig_height = 2.8
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    return fig, ax


def render_arena(arr, ax):
    """Render a 2D intensity array on the axes using the arena colormap."""
    ax.imshow(arr, cmap=ARENA_CMAP, vmin=0, vmax=15,
              interpolation='nearest', aspect='equal')
    ax.set_xlim(-0.5, ARENA_W - 0.5)
    ax.set_ylim(ARENA_H - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)


def draw_grid(ax, row_positions, col_positions, flash_size,
              edgecolor=GRID_COLOR, linewidth=GRID_LW):
    """Draw grid of rectangle outlines for flash positions."""
    for r in row_positions:
        for c in col_positions:
            rect = mpatches.Rectangle(
                (c - 1 - 0.5, r - 1 - 0.5),
                flash_size, flash_size,
                linewidth=linewidth,
                edgecolor=edgecolor,
                facecolor='none'
            )
            ax.add_patch(rect)


def draw_rect_outline(ax, row_start, row_end, col_start, col_end,
                      edgecolor=OUTLINE_COLOR, linewidth=OUTLINE_LW):
    """Draw a rectangle outline using 1-indexed inclusive coordinates."""
    width = col_end - col_start + 1
    height = row_end - row_start + 1
    rect = mpatches.Rectangle(
        (col_start - 1 - 0.5, row_start - 1 - 0.5),
        width, height,
        linewidth=linewidth,
        edgecolor=edgecolor,
        facecolor='none'
    )
    ax.add_patch(rect)


def save_figure(fig, name):
    """Save figure as both SVG and PDF."""
    svg_path = os.path.join(OUTPUT_DIR, f"{name}.svg")
    pdf_path = os.path.join(OUTPUT_DIR, f"{name}.pdf")
    fig.savefig(svg_path, format='svg', bbox_inches='tight', transparent=True)
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f"  Saved: {name}.svg, {name}.pdf")


# ============================================================================
# Protocol 2 helper: compute the 30x30 area bounds
# ============================================================================

def get_p2_area():
    """Get the 30x30 Protocol 2 area bounds for the example RF center."""
    return centered_square(RF_CENTER_X, RF_CENTER_Y, P2_CROP_SIZE)


# ============================================================================
# Figure 1: Protocol 1 - 12x12 pixel flash grid
# ============================================================================

def figure1_12px_flash_grid():
    """Protocol 1: Full arena with 12x12 non-overlapping flash grid on LHS."""
    arr = create_arena_array(OFF_INTENSITY)
    fill_region(arr, P1_ROW_START, P1_ROW_END, P1_COL_START, P1_COL_END,
                BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    flash_size = 12
    row_pos = compute_flash_positions(P1_ROW_START, P1_ROW_END, flash_size, 0)
    col_pos = compute_flash_positions(P1_COL_START, P1_COL_END, flash_size, 0)
    draw_grid(ax, row_pos, col_pos, flash_size)

    save_figure(fig, 'fig1_protocol1_12px_flash_grid')


# ============================================================================
# Figure 2: Protocol 1 - 6x6 pixel flash grid
# ============================================================================

def figure2_6px_flash_grid():
    """Protocol 1: Full arena with 6x6 non-overlapping flash grid on LHS."""
    arr = create_arena_array(OFF_INTENSITY)
    fill_region(arr, P1_ROW_START, P1_ROW_END, P1_COL_START, P1_COL_END,
                BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    flash_size = 6
    row_pos = compute_flash_positions(P1_ROW_START, P1_ROW_END, flash_size, 0)
    col_pos = compute_flash_positions(P1_COL_START, P1_COL_END, flash_size, 0)
    draw_grid(ax, row_pos, col_pos, flash_size)

    save_figure(fig, 'fig2_protocol1_6px_flash_grid')


# ============================================================================
# Figure 3: Protocol 2 - Example 6x6 flash with 30x30 area outline
# ============================================================================

def figure3_example_flash_30px():
    """Protocol 2: 6x6 grid with one example flash filled and 30x30 outline."""
    arr = create_arena_array(OFF_INTENSITY)
    fill_region(arr, P1_ROW_START, P1_ROW_END, P1_COL_START, P1_COL_END,
                BKG_INTENSITY)

    # Fill one example 6x6 flash with dark (black)
    # Choose flash at row=19, col=65 (near center of LHS)
    example_row = 19
    example_col = 65
    fill_region(arr, example_row, example_row + 5, example_col, example_col + 5,
                OFF_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    # Draw 6x6 grid
    flash_size = 6
    row_pos = compute_flash_positions(P1_ROW_START, P1_ROW_END, flash_size, 0)
    col_pos = compute_flash_positions(P1_COL_START, P1_COL_END, flash_size, 0)
    draw_grid(ax, row_pos, col_pos, flash_size)

    # Draw 30x30 area outline centered on the example flash
    r1, r2, c1, c2 = get_p2_area()
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    save_figure(fig, 'fig3_protocol2_example_flash_30px')


# ============================================================================
# Figure 4: Protocol 2 - 6x6 flash grid with 50% overlap in 30x30 area
# ============================================================================

def figure4_6px_50pct_overlap():
    """Protocol 2: 6x6 flash grid with 50% overlap in 30x30 area."""
    arr = create_arena_array(BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    r1, r2, c1, c2 = get_p2_area()
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    flash_size = P2_FLASH_6PX
    row_pos = compute_flash_positions(r1, r2, flash_size, P2_OVERLAP)
    col_pos = compute_flash_positions(c1, c2, flash_size, P2_OVERLAP)
    draw_grid(ax, row_pos, col_pos, flash_size)

    save_figure(fig, 'fig4_protocol2_6px_50pct_overlap')


# ============================================================================
# Figure 5: Protocol 2 - 4x4 flash grid with 50% overlap in 30x30 area
# ============================================================================

def figure5_4px_50pct_overlap():
    """Protocol 2: 4x4 flash grid with 50% overlap in 30x30 area."""
    arr = create_arena_array(BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    r1, r2, c1, c2 = get_p2_area()
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    flash_size = P2_FLASH_4PX
    row_pos = compute_flash_positions(r1, r2, flash_size, P2_OVERLAP)
    col_pos = compute_flash_positions(c1, c2, flash_size, P2_OVERLAP)
    draw_grid(ax, row_pos, col_pos, flash_size)

    save_figure(fig, 'fig5_protocol2_4px_50pct_overlap')


# ============================================================================
# Figure 6: Protocol 2 - Single vertical bar with sweep arrow
# ============================================================================

def figure6_vertical_bar_sweep():
    """Protocol 2: Vertical dark bar centered in 30x30 area with rightward arrow."""
    arr = create_arena_array(BKG_INTENSITY)

    r1, r2, c1, c2 = get_p2_area()

    # Center the 4px bar horizontally within the 30x30 area
    area_center_col = (c1 + c2) / 2.0
    bar_col_start = int(area_center_col - BAR_WIDTH / 2 + 1)
    bar_col_end = bar_col_start + BAR_WIDTH - 1

    fill_region(arr, r1, r2, bar_col_start, bar_col_end, OFF_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    # Arrow pointing right from bar
    arrow_y = (r1 + r2) / 2.0 - 1  # 0-indexed center row
    arrow_x_start = bar_col_end - 1 + 1.0  # just right of bar (0-indexed)
    arrow_x_end = arrow_x_start + 8.0
    ax.annotate('', xy=(arrow_x_end, arrow_y),
                xytext=(arrow_x_start, arrow_y),
                arrowprops=dict(arrowstyle='->', color='lightgrey', lw=1.5))

    save_figure(fig, 'fig6_protocol2_vertical_bar_sweep')


# ============================================================================
# Figure 7: Protocol 2 - Single diagonal bar
# ============================================================================

def figure7_diagonal_bar():
    """Protocol 2: Dark bar rotated ~45 degrees spanning 30x30 diagonal."""
    arr = create_arena_array(BKG_INTENSITY)

    r1, r2, c1, c2 = get_p2_area()

    # Draw diagonal bar by setting pixels within the 30x30 area
    # Bar is 4px wide, rotated 45 degrees, spanning the diagonal
    center_row = (r1 + r2) / 2.0
    center_col = (c1 + c2) / 2.0
    half_w = BAR_WIDTH / 2.0

    for r in range(r1, r2 + 1):
        for c in range(c1, c2 + 1):
            dr = r - center_row
            dc = c - center_col
            # Perpendicular distance to the 45-degree line (y = x) through center
            dist = abs(dr - dc) / np.sqrt(2)
            if dist < half_w:
                arr[r - 1, c - 1] = OFF_INTENSITY

    fig, ax = create_figure()
    render_arena(arr, ax)

    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    save_figure(fig, 'fig7_protocol2_diagonal_bar')


# ============================================================================
# Figure 8: Protocol 2 - 11 bar flash positions
# ============================================================================

def figure8_bar_flash_positions():
    """Protocol 2: Center bar filled dark, 10 flanking bar outlines."""
    arr = create_arena_array(BKG_INTENSITY)

    r1, r2, c1, c2 = get_p2_area()

    # Center bar column (same as figure 6)
    area_center_col = (c1 + c2) / 2.0
    center_col_start = int(area_center_col - BAR_WIDTH / 2 + 1)

    # Fill center bar dark
    fill_region(arr, r1, r2, center_col_start, center_col_start + BAR_WIDTH - 1,
                OFF_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    # Draw 11 bar positions: center + 5 left + 5 right, shifted by 2px each
    bar_height = r2 - r1 + 1
    for i in range(N_BAR_POSITIONS):
        offset = (i - N_FLANK) * BAR_POSITION_STEP
        col_start = center_col_start + offset

        if i == N_FLANK:
            # Center bar is already filled in the raster; draw its outline too
            draw_rect_outline(ax, r1, r2, col_start,
                              col_start + BAR_WIDTH - 1,
                              edgecolor=GRID_COLOR, linewidth=GRID_LW)
        else:
            # Flanking bars: outline only
            draw_rect_outline(ax, r1, r2, col_start,
                              col_start + BAR_WIDTH - 1,
                              edgecolor=GRID_COLOR, linewidth=GRID_LW)

    save_figure(fig, 'fig8_protocol2_bar_flash_positions')


# ============================================================================
# Main
# ============================================================================

def main():
    """Generate all 8 schematic figures."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Generating LED arena stimulus schematics...")
    print(f"Output directory: {OUTPUT_DIR}\n")

    print("Protocol 1:")
    figure1_12px_flash_grid()
    figure2_6px_flash_grid()

    print("\nProtocol 2:")
    figure3_example_flash_30px()
    figure4_6px_50pct_overlap()
    figure5_4px_50pct_overlap()
    figure6_vertical_bar_sweep()
    figure7_diagonal_bar()
    figure8_bar_flash_positions()

    print(f"\nAll 8 figures saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
