function [bf_col, positions_order] = map_PD_to_bar_flash_column(pd_angle_rad)
    % Map a preferred direction angle from bar sweep analysis to the
    % corresponding bar flash orientation column and determine leading-to-
    % trailing position ordering.
    %
    % Inputs
    % ------
    %   pd_angle_rad : float
    %       Preferred direction angle in radians (0 to 2*pi) from bar sweep
    %       analysis.
    %
    % Outputs
    % -------
    %   bf_col : int
    %       Column index (1-8) in the bar flash data corresponding to the
    %       orientation whose position axis aligns with the PD.
    %
    %   positions_order : array [1x11]
    %       Indices 1:11 ordered from leading to trailing position along the
    %       PD axis.
    %
    % Bar flash orientation mapping (from generate_bar_flash_stimulus_xy.m):
    %   Column 1: pi     (180 deg)
    %   Column 2: pi/8   (22.5 deg)
    %   Column 3: 2pi/8  (45 deg)
    %   Column 4: 3pi/8  (67.5 deg)
    %   Column 5: 4pi/8  (90 deg)
    %   Column 6: 5pi/8  (112.5 deg)
    %   Column 7: 6pi/8  (135 deg)
    %   Column 8: 7pi/8  (157.5 deg)

    % Bar flash column orientations (the angle at which the bar is oriented)
    bf_orientations = [pi, pi/8, 2*pi/8, 3*pi/8, 4*pi/8, 5*pi/8, 6*pi/8, 7*pi/8];

    % A bar sweep moving at pd_angle_rad has bars oriented perpendicular to
    % the motion direction.
    bar_orient = mod(pd_angle_rad - pi/2, pi); % normalise to [0, pi)

    % Normalise bf_orientations to [0, pi) for comparison
    bf_orient_mod = mod(bf_orientations, pi);

    % Find the closest bar flash orientation
    diffs = abs(bf_orient_mod - bar_orient);
    diffs = min(diffs, pi - diffs); % handle wrap-around at pi
    [~, bf_col] = min(diffs);

    % Determine leading-to-trailing position ordering.
    % The position axis is perpendicular to the bar orientation.
    pos_axis_angle = bf_orientations(bf_col) + pi/2;

    % Check if PD aligns with the positive direction of the position axis
    cos_alignment = cos(pd_angle_rad - pos_axis_angle);

    if cos_alignment >= 0
        positions_order = 1:11;   % position 1 is leading
    else
        positions_order = 11:-1:1; % position 11 is leading
    end

    fprintf('  PD angle: %.1f deg -> bar flash column %d (bar orient %.1f deg), positions %d->%d\n', ...
        rad2deg(pd_angle_rad), bf_col, rad2deg(bf_orientations(bf_col)), ...
        positions_order(1), positions_order(end));

end
