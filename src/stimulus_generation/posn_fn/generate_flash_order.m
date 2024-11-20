function outputSequence = generate_flash_order(screen_size, pixel_size)

% Generate the order in which flashes should be presented in the RF
% stimulus. The screen is split up into a grid based on the size of the
% flashes. This grid is then split into quadrants and flashes are shown
% sequenctially in the same position within each quadrant in order.

% This works with a pattern that has the flashes ordered going down the
% ROWS first, then across the COLUMNS. 

% Required inputs:
% - pixel_size 
% - size of screen used (in pixels)
% - additionally in future, if the grid is overlapping or non-overlapping.

% ______________________________________________________________

    % Determine the size of the grid of flashes
    grid_cols = screen_size(2)/pixel_size; % Total width of the image
    grid_rows = screen_size(1)/pixel_size; % Total height of the image
    
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














