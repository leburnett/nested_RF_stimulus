function cmap = create_red_blue_cmap()

    % Build custom 11×3 colormap
    cmap = zeros(11,3);
    
    % 1–5: red → white
    startRed = [1 0 0];  % pure red
    endWhite = [1 1 1];  % white
    for k = 1:6
        t = (k-1)/(6-1);           % interpolate 0→1 across positions 1–5
        cmap(k,:) = (1-t)*startRed + t*endWhite;
    end
    
    % 7–11: white → dark blue
    lightBlue = [0.6 0.8 1]; 
    endBlue = [0 0 0.4];           % dark blue
    for k = 7:11
        t = (k-4)/(11-4);          % interpolate 0→1 across positions 6–11
        cmap(k,:) = (1-t)*lightBlue + t*endBlue;
    end

end 