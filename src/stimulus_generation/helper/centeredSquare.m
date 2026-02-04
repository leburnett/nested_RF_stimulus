function [row_start, row_end, col_start, col_end] = centeredSquare(x, y, square_size)
% CENTEREDSQUARE  Calculate screen region boundaries centered on a point.
%
%   [ROW_START, ROW_END, COL_START, COL_END] = CENTEREDSQUARE(X, Y, SQUARE_SIZE)
%   computes the pixel boundaries of a square region centered on [x,y],
%   ensuring the region stays within the valid screen area.
%
%   INPUTS:
%     x           - Horizontal pixel coordinate (column, 17-192)
%     y           - Vertical pixel coordinate (row, 1-48)
%     square_size - Width/height of the square region in pixels
%
%   OUTPUTS:
%     row_start - First row (y) of the region (1-48)
%     row_end   - Last row (y) of the region (1-48)
%     col_start - First column (x) of the region (17-192)
%     col_end   - Last column (x) of the region (17-192)
%
%   SCREEN DIMENSIONS:
%     The G4 arena display is 192x48 pixels (12 panels x 3 panels).
%     Valid display area: columns 17-192, rows 1-48.
%     Column 1-16 is typically not used (panel configuration).
%
%   EDGE HANDLING:
%     If centering the square would extend beyond screen boundaries,
%     the region is shifted to fit within the valid screen area.
%     This ensures stimuli remain fully visible even for RF centers
%     near screen edges.
%
%   EXAMPLE:
%     % Get 30x30 pixel region centered at position [100, 24]
%     [r1, r2, c1, c2] = centeredSquare(100, 24, 30);
%     % Returns: r1=9, r2=38, c1=85, c2=114
%
%   See also GENERATE_FLASH_STIMULUS_XY, GENERATE_BAR_STIMULUS_XY
%___________________________________________________________________________
    
    % Screen dimensions
    screen_width_start = 17;
    screen_width_end = 192;
    screen_height_start = 1;
    screen_height_end = 48;

    % Square size and half-size for centering
    % square_size = 30;
    half_size = square_size / 2;

    % Calculate row and column start and end positions
    row_start = max(screen_height_start, min(y - half_size + 1, screen_height_end - square_size + 1));
    row_end = row_start + square_size - 1;

    col_start = max(screen_width_start, min(x - half_size + 1, screen_width_end - square_size + 1));
    col_end = col_start + square_size - 1;

    % Output the square boundaries
    % fprintf('Square bounds: [%d, %d, %d, %d]\n', row_start, row_end, col_start, col_end);
end