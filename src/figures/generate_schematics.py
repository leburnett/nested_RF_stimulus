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
from matplotlib.collections import PatchCollection
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

# Grid line style — use light grey for visibility on green/black background
GRID_COLOR = (0.7, 0.7, 0.7)
GRID_LW = 0.8
GRID_LW_THIN = 0.5
OUTLINE_COLOR = (0.8, 0.8, 0.8)
OUTLINE_LW = 1.5

# LED pixel circle style
PIXEL_RADIUS = 0.35
PIXEL_EDGE_LW = 0.15

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
    fig_height = fig_width * ARENA_H / ARENA_W
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    return fig, ax


def render_arena(arr, ax):
    """Render arena as a green/black raster with LED pixel circles on top.

    First draws the intensity array as a raster image using imshow,
    then overlays small circles at every pixel position to show the
    LED pixel structure.
    """
    # Draw the raster background
    ax.imshow(arr, cmap=ARENA_CMAP, vmin=0, vmax=15, aspect='equal',
              interpolation='nearest', origin='upper')

    # Draw LED pixel circles on top
    circles = []
    for r in range(ARENA_H):
        for c in range(ARENA_W):
            circles.append(mpatches.Circle((c, r), PIXEL_RADIUS))

    collection = PatchCollection(
        circles,
        facecolors='none',
        edgecolors=[(0, 0.15, 0)] * len(circles),
        linewidths=PIXEL_EDGE_LW,
    )
    ax.add_collection(collection)

    ax.set_xlim(-0.5, ARENA_W - 0.5)
    ax.set_ylim(ARENA_H - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.set_aspect('equal')


def draw_grid(ax, row_positions, col_positions, flash_size,
              edgecolor=GRID_COLOR, linewidth=GRID_LW, linestyle='-'):
    """Draw grid of rectangle outlines for flash positions."""
    for r in row_positions:
        for c in col_positions:
            rect = mpatches.Rectangle(
                (c - 1 - 0.5, r - 1 - 0.5),
                flash_size, flash_size,
                linewidth=linewidth,
                edgecolor=edgecolor,
                facecolor='none',
                linestyle=linestyle
            )
            ax.add_patch(rect)


def draw_rect_outline(ax, row_start, row_end, col_start, col_end,
                      edgecolor=OUTLINE_COLOR, linewidth=OUTLINE_LW,
                      linestyle='-'):
    """Draw a rectangle outline using 1-indexed inclusive coordinates."""
    width = col_end - col_start + 1
    height = row_end - row_start + 1
    rect = mpatches.Rectangle(
        (col_start - 1 - 0.5, row_start - 1 - 0.5),
        width, height,
        linewidth=linewidth,
        edgecolor=edgecolor,
        facecolor='none',
        linestyle=linestyle
    )
    ax.add_patch(rect)


def save_figure(fig, name):
    """Save figure as both SVG and PDF."""
    svg_path = os.path.join(OUTPUT_DIR, f"{name}.svg")
    pdf_path = os.path.join(OUTPUT_DIR, f"{name}.pdf")
    fig.savefig(svg_path, format='svg', bbox_inches='tight',
                facecolor='black', edgecolor='none')
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight',
                facecolor='black', edgecolor='none')
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
    arr = create_arena_array(BKG_INTENSITY)

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
    arr = create_arena_array(BKG_INTENSITY)

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
    """Protocol 2: 6x6 grid with one example flash filled and 30x30 outline.

    The 30x30 area is centered on the chosen example flash.
    """
    arr = create_arena_array(BKG_INTENSITY)

    # Choose a 6x6 flash near center of LHS that falls on the P1 grid.
    # P1 6px grid row positions: [1, 7, 13, 19, 25, 31, 37, 43]
    # P1 6px grid col positions: [17, 23, 29, 35, 41, 47, 53, 59, 65, 71, ...]
    # Pick flash at row=19, col=65 → covers rows 19-24, cols 65-70
    # Flash center: row ~21.5, col ~67.5
    example_row = 19
    example_col = 65
    flash_center_x = example_col + 2  # col 67 (center of 65-70)
    flash_center_y = example_row + 2  # row 21 (center of 19-24)

    # Fill the example flash with dark (black)
    fill_region(arr, example_row, example_row + 5, example_col, example_col + 5,
                OFF_INTENSITY)

    # Compute 30x30 area bounds first so we can split grid drawing
    r1, r2, c1, c2 = centered_square(flash_center_x, flash_center_y,
                                      P2_CROP_SIZE)

    fig, ax = create_figure()
    render_arena(arr, ax)

    # Draw 6x6 grid over the LHS region, split into inside/outside 30x30
    flash_size = 6
    row_pos = compute_flash_positions(P1_ROW_START, P1_ROW_END, flash_size, 0)
    col_pos = compute_flash_positions(P1_COL_START, P1_COL_END, flash_size, 0)

    # A flash cell at (r, c) with size 6 is inside the 30x30 area if its
    # entire extent overlaps (even partially) with [r1..r2, c1..c2]
    outer_grid_color = (0.55, 0.55, 0.55)
    for r in row_pos:
        for c in col_pos:
            cell_r_end = r + flash_size - 1
            cell_c_end = c + flash_size - 1
            inside = (cell_r_end >= r1 and r <= r2 and
                      cell_c_end >= c1 and c <= c2)
            color = GRID_COLOR if inside else outer_grid_color
            rect = mpatches.Rectangle(
                (c - 1 - 0.5, r - 1 - 0.5),
                flash_size, flash_size,
                linewidth=GRID_LW_THIN,
                edgecolor=color,
                facecolor='none',
                linestyle='-'
            )
            ax.add_patch(rect)

    # Draw 30x30 area outline
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    save_figure(fig, 'fig3_protocol2_example_flash_30px')


# ============================================================================
# Figure 4: Protocol 2 - 6x6 flash grid with 50% overlap in 30x30 area
# ============================================================================

def figure4_6px_50pct_overlap():
    """Protocol 2: Two offset non-overlapping 6x6 grids showing 50% overlap."""
    arr = create_arena_array(BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    r1, r2, c1, c2 = get_p2_area()
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    flash_size = P2_FLASH_6PX
    step = int(flash_size * (1 - P2_OVERLAP))  # = 3

    # Sub-grid A: non-overlapping grid at full flash_size spacing (solid)
    row_pos_a = list(range(r1, r2 - flash_size + 2, flash_size))
    col_pos_a = list(range(c1, c2 - flash_size + 2, flash_size))
    draw_grid(ax, row_pos_a, col_pos_a, flash_size,
              linewidth=GRID_LW_THIN, linestyle='-')

    # Sub-grid B: same grid but offset by half a flash (dashed)
    row_pos_b = list(range(r1 + step, r2 - flash_size + 2, flash_size))
    col_pos_b = list(range(c1 + step, c2 - flash_size + 2, flash_size))
    draw_grid(ax, row_pos_b, col_pos_b, flash_size,
              linewidth=GRID_LW_THIN, linestyle='--')

    save_figure(fig, 'fig4_protocol2_6px_50pct_overlap')


# ============================================================================
# Figure 5: Protocol 2 - 4x4 flash grid with 50% overlap in 30x30 area
# ============================================================================

def figure5_4px_50pct_overlap():
    """Protocol 2: Two offset non-overlapping 4x4 grids showing 50% overlap."""
    arr = create_arena_array(BKG_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    r1, r2, c1, c2 = get_p2_area()
    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    flash_size = P2_FLASH_4PX
    step = int(flash_size * (1 - P2_OVERLAP))  # = 2

    # Sub-grid A: non-overlapping grid at full flash_size spacing (solid)
    row_pos_a = list(range(r1, r2 - flash_size + 2, flash_size))
    col_pos_a = list(range(c1, c2 - flash_size + 2, flash_size))
    draw_grid(ax, row_pos_a, col_pos_a, flash_size,
              linewidth=GRID_LW_THIN, linestyle='-')

    # Sub-grid B: same grid but offset by half a flash (dashed)
    row_pos_b = list(range(r1 + step, r2 - flash_size + 2, flash_size))
    col_pos_b = list(range(c1 + step, c2 - flash_size + 2, flash_size))
    draw_grid(ax, row_pos_b, col_pos_b, flash_size,
              linewidth=GRID_LW_THIN, linestyle='--')

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
                arrowprops=dict(arrowstyle='->', color=OUTLINE_COLOR, lw=1.5))

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
    """Protocol 2: Center bar filled, edge lines for all positions, brackets.

    Only the central bar is filled dark. For each of the 11 positions,
    left and right boundary vertical lines are drawn (alternating solid/dotted).
    Brackets showing x-extent of each bar are drawn INSIDE the arena:
    left-of-center bars above the 30x30 area (between arena top and 30x30 top),
    right-of-center + center bars below the 30x30 area (between 30x30 bottom
    and arena bottom).

    The 30x30 area spans rows 10-39 (0-indexed 9-38), leaving 9 pixel rows
    above (rows 0-8) and 9 pixel rows below (rows 39-47) for brackets.
    """
    arr = create_arena_array(BKG_INTENSITY)

    r1, r2, c1, c2 = get_p2_area()

    # Center bar column (same as figure 6)
    area_center_col = (c1 + c2) / 2.0
    center_col_start = int(area_center_col - BAR_WIDTH / 2 + 1)

    # Fill only the center bar dark in the raster
    fill_region(arr, r1, r2, center_col_start, center_col_start + BAR_WIDTH - 1,
                OFF_INTENSITY)

    fig, ax = create_figure()
    render_arena(arr, ax)

    draw_rect_outline(ax, r1, r2, c1, c2, linewidth=1.5)

    bracket_color = (0.8, 0.8, 0.8)
    bracket_lw = 0.8
    edge_line_lw = 0.7

    # Bracket layout:
    # 30x30 area: r1=10, r2=39 (1-indexed) → 0-indexed rows 9..38
    # Above 30x30 area: 0-indexed rows 0..8 (9 rows of space for 5 brackets)
    # Below 30x30 area: 0-indexed rows 39..47 (9 rows of space for 6 brackets)
    # Bracket spacing: ~1.5 px per level, tick height ~1.0 px
    bracket_spacing = 1.5
    bracket_tick_h = 1.0

    # The 30x30 area top edge in 0-indexed coords
    area_top_0 = r1 - 1 - 0.5    # top edge of 30x30 area
    area_bot_0 = r2 - 0.5        # bottom edge of 30x30 area

    for i in range(N_BAR_POSITIONS):
        offset = (i - N_FLANK) * BAR_POSITION_STEP
        col_start = center_col_start + offset
        col_end = col_start + BAR_WIDTH - 1

        # Alternate solid/dotted for edge lines
        if i % 2 == 0:
            ls = '-'
        else:
            ls = ':'

        # Draw left and right boundary vertical lines (full height of 30x30)
        bx_left_0 = col_start - 1 - 0.5   # 0-indexed, left pixel edge
        bx_right_0 = col_end - 0.5         # 0-indexed, right pixel edge

        ax.plot([bx_left_0, bx_left_0], [area_top_0, area_bot_0],
                color=GRID_COLOR, linewidth=edge_line_lw, linestyle=ls)
        ax.plot([bx_right_0, bx_right_0], [area_top_0, area_bot_0],
                color=GRID_COLOR, linewidth=edge_line_lw, linestyle=ls)

        # Draw brackets inside the arena
        dist_from_center = abs(i - N_FLANK)

        if i < N_FLANK:
            # Left-of-center bars: brackets ABOVE the 30x30 area
            # Stack from the 30x30 top edge upward into arena rows 0-8
            # Level 1 (closest to center) nearest the 30x30 edge, level 5 furthest
            level = dist_from_center  # 1..5
            bracket_y_bottom = area_top_0 - 0.5 - (level - 1) * bracket_spacing
            bracket_y_top = bracket_y_bottom - bracket_tick_h
            # Inverted U bracket above the 30x30 area
            ax.plot([bx_left_0, bx_left_0],
                    [bracket_y_bottom, bracket_y_top],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')
            ax.plot([bx_right_0, bx_right_0],
                    [bracket_y_bottom, bracket_y_top],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')
            ax.plot([bx_left_0, bx_right_0],
                    [bracket_y_top, bracket_y_top],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')
        else:
            # Center and right-of-center bars: brackets BELOW the 30x30 area
            # Stack from the 30x30 bottom edge downward into arena rows 39-47
            if i == N_FLANK:
                level = 0  # center is closest
            else:
                level = dist_from_center  # 1..5
            bracket_y_top = area_bot_0 + 0.5 + level * bracket_spacing
            bracket_y_bottom = bracket_y_top + bracket_tick_h
            # U bracket below the 30x30 area
            ax.plot([bx_left_0, bx_left_0],
                    [bracket_y_top, bracket_y_bottom],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')
            ax.plot([bx_right_0, bx_right_0],
                    [bracket_y_top, bracket_y_bottom],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')
            ax.plot([bx_left_0, bx_right_0],
                    [bracket_y_bottom, bracket_y_bottom],
                    color=bracket_color, linewidth=bracket_lw,
                    solid_capstyle='butt')

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
