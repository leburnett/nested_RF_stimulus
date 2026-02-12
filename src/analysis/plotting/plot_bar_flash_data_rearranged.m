function f = plot_bar_flash_data_rearranged(meanData)

    % Inputs:
    % meanData : cell array containing vectors, [nPos x nOrient]

    % Plots the mean responses to bar flashes. Rearranged so that the
    % position with the strongest response is positioned in the middle
    % (6th) column. The background colour corresponds to the bar position.

    % create_red_blue_cmap : function that returns 11x3 colormap
    cmap = create_red_blue_cmap();   % 11×3, index 1..11 correspond to bar positions
    nPos = 11;                    % positions
    nOrient = 8;               % orientations
    
    % --- Precompute an overall Y max for consistent y-limits ---
    max_overall = -inf;
    for o = 1:nOrient
        for p = 1:nPos
            d = meanData{p, o};
            if ~isempty(d)
                max_overall = max(max_overall, max(d(:)));
            end
        end
    end
    if ~isfinite(max_overall)
        error('meanData appears empty.');
    end
    
    % --- Plot settings ---
    y_base   = -70;    % baseline for background rectangle
    y_height = 40;     % height of background rectangle
    line_w   = 1.75;      % linewidth
    txt_sz   = 10;      % title font size
    
    % --- Build figure with 1 row per orientation, 11 columns (positions) ---
    figure;
    tiledlayout(nOrient, 11, 'TileSpacing','compact','Padding','compact');
    
    for o = 1:nOrient
        % 1) Find the position with the largest response in this orientation
        peakVals = zeros(1, nPos);
        for p = 1:nPos
            d = meanData{p, o};
            peakVals(p) = max(d);     % change to max(abs(d)) if that's your metric
        end
        [~, best_pos] = max(peakVals);
    
        % 2) Make the column order so best_pos lands in column 6, others retain order
        order = 1:nPos;                         % original left→right positions
        rot   = mod(6 - best_pos, nPos);        % rotation to place best at col 6
        colOrder = circshift(order, rot);       % what we will PLOT left→right
    
        % 3) Plot this orientation’s row, following colOrder
        for c = 1:nPos
            pos = colOrder(c);                  % original position being plotted in column c
            d   = meanData{pos, o};
    
            nexttile((o-1)*11 + c);
            % Background color by ORIGINAL position (before reordering)
            rectangle('Position', [0, y_base, numel(d), y_height], ...
                      'FaceColor', cmap(pos,:), 'EdgeColor','none', 'FaceAlpha', 0.6); hold on;
    
            plot(d, 'k', 'LineWidth', line_w);
            xlim([0, numel(d)]);
            ylim([y_base, max_overall*0.9]);
            title(sprintf('Pos %d', pos), 'FontSize', txt_sz);
            axis tight; box off;
    
            if o == 1 && c ==1 
                ylabel('Voltage (mV)')
                xticks([0 5000 10000])
                xticklabels({'0', '0.5', '1'})
                xlabel('Time (s)')
            else
                ax = gca;
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
            end 
        end
    end
    f = gcf;
    f.Position = [10 84 1773 961];

end 




