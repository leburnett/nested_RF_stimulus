function outputSequence = generate_flash_order(fl_rows, fl_cols)
% GENERATE_FLASH_ORDER  Create quadrant-interleaved flash presentation order.
%
%   OUTPUTSEQUENCE = GENERATE_FLASH_ORDER(FL_ROWS, FL_COLS)
%   generates a sequence of frame indices that presents flashes in a
%   spatially distributed order to minimize adaptation effects between
%   adjacent flashes.
%
%   INPUTS:
%     fl_rows - Number of rows in the flash grid
%     fl_cols - Number of columns in the flash grid
%
%   OUTPUT:
%     outputSequence - 1 x (fl_rows*fl_cols) array of frame indices
%                      specifying the order to present flashes
%
%   ALGORITHM:
%     The flash grid is divided into four quadrants:
%       Q1 (top-left)     Q2 (top-right)
%       Q3 (bottom-left)  Q4 (bottom-right)
%
%     For each position within a quadrant, the function cycles through
%     all four quadrants before moving to the next position. This ensures
%     consecutive flashes are spatially separated, reducing local
%     adaptation artifacts in the neural response.
%
%   PATTERN COMPATIBILITY:
%     Designed for patterns with column-major ordering (flashes ordered
%     going down rows first, then across columns). Uses sub2ind() for
%     conversion between [row,col] and linear frame indices.
%
%   EXAMPLE:
%     % For a 14x14 flash grid (196 flashes)
%     order = generate_flash_order(14, 14);
%     % order = [1, 8, 99, 106, 2, 9, 100, 107, ...]
%     % Cycles through all 4 quadrants for each position
%
%   See also GENERATE_FLASH_FUNCTION, GENERATE_FUNC_FOR_FLASH
% ______________________________________________________________

    % Determine the size of the grid of flashes
    grid_cols = fl_cols; % Total width of the image
    grid_rows = fl_rows; % Total height of the image
    
    % Determine the size of the grid quadrants.
    quadrant_cols = grid_cols / 2; % Width of each quadrant
    quadrant_rows = grid_rows / 2; % Height of each quadrant
    
    % Initialize the output sequence
    outputSequence = zeros(1, grid_cols * grid_rows);
    index = 1;
    
    % Loop through all the positions within each quadrant
    for row = 1:quadrant_rows
        for col = 1:quadrant_cols
            % Define the top-left pixels of the four quadrants
            % Quadrant 1 (top-left)
            outputSequence(index) = sub2ind([grid_rows, grid_cols], row, col);
            index = index + 1;
    
            % Quadrant 2 (top-right)
            outputSequence(index) = sub2ind([grid_rows, grid_cols], row, col + quadrant_cols);
            index = index + 1;
    
            % Quadrant 3 (bottom-left)
            outputSequence(index) = sub2ind([grid_rows, grid_cols], row + quadrant_rows, col);
            index = index + 1;
    
            % Quadrant 4 (bottom-right)
            outputSequence(index) = sub2ind([grid_rows, grid_cols], row + quadrant_rows, col + quadrant_cols);
            index = index + 1;
        end
    end

end 














