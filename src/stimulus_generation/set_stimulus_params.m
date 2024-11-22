% Generate both the pattern and position function for RF flash stimuli.   

%% Region of arena on which to present the stimulus. 

%  1 - For the 1/2 screen area. 6 panels. RHS of screen when look from the
%  front. LHS from fly's perspective. For recording from left hemisphere. 

disp_x1 = 17;
disp_x2 = 112;
disp_y1 = 1;
disp_y2 = 48;

%  2 - For the 1/2 screen area. 6 panels. LHS of screen when look from the
%  front. RHS from fly's perspective. For recording from right hemisphere. 

% disp_x1 = 101;
% disp_x2 = 196;
% disp_y1 = 1;
% disp_y2 = 48;

% 3 - 30 x 30 pixel square

% All possible 30 x 30 pixel squares in the area 17:196 , 1:48. 

% disp_x1 = 17;
% disp_x2 = 46;
% disp_y1 = 1;
% disp_y2 = 30;


%% Other parameters which don't change with the area of the screen that the
% pattern is being presented. 

params.px_intensity = [6, 0, 15]; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
params.flash_sz_px = 6; % Size of flash in pixels
params.overlap = 0; % Overlap of flashes - between 0 and 1. 

params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
params.interval_dur = 0.05; % duration of interval background screen in seconds.
params.flash_dur = 0.2; % duration of flash in seconds.

generate_stimulus(params)






% Screen dimensions and limits
screen_width_start = 17;
screen_width_end = 196;
screen_height_start = 1;
screen_height_end = 48;

% Square size
square_size = 30;

% Initialize array to store possible square positions
positions = [];

% Loop through all possible row and column start positions
for row_start = screen_height_start:(screen_height_end - square_size + 1)
    row_end = row_start + square_size - 1;  % Calculate row end
    for col_start = screen_width_start:(screen_width_end - square_size + 1)
        col_end = col_start + square_size - 1;  % Calculate column end
        
        % Append the position to the list
        positions = [positions; row_start, row_end, col_start, col_end];
    end
end

% Display all possible positions
disp(positions);


function [row_start, row_end, col_start, col_end] = centeredSquare(x, y)
    % Screen dimensions
    screen_width_start = 17;
    screen_width_end = 196;
    screen_height_start = 1;
    screen_height_end = 48;

    % Square size and half-size for centering
    square_size = 30;
    half_size = square_size / 2;

    % Calculate row and column start and end positions
    row_start = max(screen_height_start, min(y - half_size + 1, screen_height_end - square_size + 1));
    row_end = row_start + square_size - 1;

    col_start = max(screen_width_start, min(x - half_size + 1, screen_width_end - square_size + 1));
    col_end = col_start + square_size - 1;

    % Output the square boundaries
    fprintf('Square bounds: [%d, %d, %d, %d]\n', row_start, row_end, col_start, col_end);
end
