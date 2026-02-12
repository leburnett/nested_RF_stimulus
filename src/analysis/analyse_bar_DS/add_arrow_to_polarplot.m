function add_arrow_to_polarplot(magnitude, angl, col)
% ADD_ARROW_TO_POLARPLOT  Overlay a vector arrow on an existing polar plot.
%
%   ADD_ARROW_TO_POLARPLOT(MAGNITUDE, ANGL, COL) draws an arrow from the
%   origin to the specified polar coordinates, useful for indicating
%   the preferred direction from a vector sum analysis.
%
%   INPUTS:
%     magnitude - Length of the arrow (in plot units)
%     angl      - Direction of arrow in radians (0 = East, pi/2 = North)
%     col       - Color specification (RGB triplet or color name)
%
%   USAGE:
%     Must be called after creating a polarplot. The arrow will be
%     overlaid on the existing polar axes.
%
%   ARROW DESIGN:
%     - Main shaft: thick line from origin
%     - Arrowhead: triangular, 30 degree angle
%     - Arrowhead size proportional to arrow length (1/15)
%
%   EXAMPLE:
%     polarplot(theta, rho);
%     hold on
%     add_arrow_to_polarplot(30, pi/4, [0.3 0.3 0.3])  % Arrow pointing NE
%
%   See also PLOT_POLAR_WITH_ARROW, VECTOR_SUM_POLAR 

    %%%%arrow head %%%%
    arrowhead_length    = magnitude/15; % arrow head length relative to resultant_length
    num_arrowlines = 100;
    arrowhead_angle = deg2rad(30); % degrees

    %%%%arrow tip coordinates %%%%
    t1 = repmat(angl,1,num_arrowlines);
    r1 = repmat(magnitude,1,num_arrowlines);

    %%%%arrow base coordinates %%%%
    b = arrowhead_length.*tan(linspace(0,arrowhead_angle,num_arrowlines/2));
    theta = atan(b./(magnitude-arrowhead_length));
    pre_t2 = [theta, -theta];
    r2 = (magnitude-arrowhead_length)./cos(pre_t2);
    t2 = t1(1)+pre_t2;
    %%%%plot %%%%
    
    polarplot([t1(1) t1(1)],[0 r1(1)-0.9*arrowhead_length],'Color', col,'linewidth', 3)
    hold on
    polarplot([t1; t2],[r1; r2], 'Color', col)

end 