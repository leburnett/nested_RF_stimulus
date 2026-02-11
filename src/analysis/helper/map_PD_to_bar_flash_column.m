function [bf_col, positions_order] = map_PD_to_bar_flash_column(pd_data_row, pd_angle_rad)
    % Map a preferred direction from bar sweep analysis to the corresponding
    % bar flash orientation column.
    %
    % The sweep patterns (files 0003-0010) and bar flash columns (1-8) use
    % the same sequence of rotation angles applied to a vertical bar:
    %   [pi, pi/8, 2pi/8, 3pi/8, 4pi/8, 5pi/8, 6pi/8, 7pi/8]
    %
    % Therefore the mapping is direct:
    %   Pattern file index = ceil(data_row / 2) + 2
    %   bf_col = pattern_file_index - 2 = ceil(data_row / 2)
    %
    % Inputs
    % ------
    %   pd_data_row : int
    %       Data row (1-16) from plot_order for the preferred direction.
    %       Odd rows = forward, even rows = flip. Rows 1-2 share pattern 3,
    %       rows 3-4 share pattern 4, etc.
    %
    %   pd_angle_rad : float
    %       Preferred direction angle in radians (0 to 2*pi). Used only for
    %       determining leading-to-trailing position ordering.
    %
    % Outputs
    % -------
    %   bf_col : int
    %       Column index (1-8) in the bar flash data.
    %
    %   positions_order : array [1x11]
    %       Indices 1:11 ordered from leading to trailing position along the
    %       PD axis.

    % Direct mapping: data rows 1-2 -> bf_col 1, rows 3-4 -> bf_col 2, etc.
    bf_col = ceil(pd_data_row / 2);

    % Rotation angles used to generate both sweep and flash patterns
    rotation_angles = [pi, pi/8, 2*pi/8, 3*pi/8, 4*pi/8, 5*pi/8, 6*pi/8, 7*pi/8];

    % Determine leading-to-trailing position ordering.
    % The bar flash positions are arranged along the axis perpendicular to
    % the bar orientation. We need to figure out whether position 1 or
    % position 11 is the "leading" position (closest to where the bar
    % enters from in the PD).
    %
    % For now, use the forward/flip status: forward (odd data rows) and
    % flip (even data rows) sweep in opposite directions. The forward
    % direction plays frames 11->62, flip plays 62->11.
    is_flip = mod(pd_data_row, 2) == 0;

    if ~is_flip
        positions_order = 1:11;   % forward: position 1 is leading
    else
        positions_order = 11:-1:1; % flip: position 11 is leading
    end

    fprintf('  PD data row: %d -> bf_col %d (rotation angle %.1f deg), positions %d->%d\n', ...
        pd_data_row, bf_col, rad2deg(rotation_angles(bf_col)), ...
        positions_order(1), positions_order(end));

end
