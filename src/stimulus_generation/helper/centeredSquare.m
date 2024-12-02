function [row_start, row_end, col_start, col_end] = centeredSquare(x, y, square_size)
% Find the coordinates of the subsection of the screen within which the
% stimulus should be presented, centred on the coordinates [x,y].

% [x,y] are coordinates referring to the pixels of the screen. This coordinate 
% is found by working backwards from what flash stimulus evokes the greatest
% response in protocol 1. 
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