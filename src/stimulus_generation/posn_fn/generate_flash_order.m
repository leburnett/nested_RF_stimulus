function outputSequence = generate_flash_order(screen_size, pixel_size)

% Generate the order in which flashes should be presented in the RF
% stimulus. This works with a pattern that has the flashes ordered going down the
% ROWS first, then across the COLUMNS. 

% Required inputs:
% - pixel_size 
% - size of screen used (in pixels)
% - additionally in future, if the grid is overlapping or non-overlapping.

% ______________________________________________________________

% Define the dimensions of the image and quadrants
imageWidth = screen_size(2)/pixel_size; % Total width of the image
imageHeight = screen_size(1)/pixel_size; % Total height of the image
quadrantWidth = imageWidth / 2; % Width of each quadrant
quadrantHeight = imageHeight / 2; % Height of each quadrant

% Initialize the output sequence
outputSequence = zeros(1, imageWidth * imageHeight);
index = 1;

% Loop through all the positions within each quadrant
for row = 1:quadrantHeight
    for col = 1:quadrantWidth
        % Define the top-left pixels of the four quadrants
        % Quadrant 1 (top-left)
        outputSequence(index) = sub2ind([imageHeight, imageWidth], row, col);
        index = index + 1;

        % Quadrant 2 (top-right)
        outputSequence(index) = sub2ind([imageHeight, imageWidth], row, col + quadrantWidth);
        index = index + 1;

        % Quadrant 3 (bottom-left)
        outputSequence(index) = sub2ind([imageHeight, imageWidth], row + quadrantHeight, col);
        index = index + 1;

        % Quadrant 4 (bottom-right)
        outputSequence(index) = sub2ind([imageHeight, imageWidth], row + quadrantHeight, col + quadrantWidth);
        index = index + 1;
    end
end















