function add_arrow_to_polarplot(magnitude, angl, col)
% Plot an arrow over the top of a polar plot, with the arrow length 'magnitude' and direction 'angl'. 

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