
% Analysing responses to bar stimuli

f_data = Log.ADC.Volts(1, :);
v_data = Log.ADC.Volts(2, :);

median_voltage = median(v_data)*10;

% Manually finding the timepoints at which the stimuli occur. 

% Rep 1
rep1_rng = [1959080, 5329750];
rep2_rng = [7287720, 10656500];
rep3_rng = [12615800, 15984550]; 

% 16 directions within each rep. 

%% Rep 1 
st_val = rep1_rng(1);
end_val = rep1_rng(2);

frames_rep1 = f_data(st_val:end_val);

r1_min = islocalmin(frames_rep1);
r1_min_vals = find(r1_min==1);
r1_st = [st_val, r1_min_vals+st_val, end_val];

r1_max = islocalmax(frames_rep1);
r1_max_vals = find(r1_max==1);
r1_nd = r1_max_vals+st_val;

all_idxs = [r1_st, r1_nd];
idxs = sort(all_idxs);

% TEST
% figure; plot(a)
% hold on
% for i = 1:49
% plot([idxs(i), idxs(i)], [0 100], 'r')
% end 

%% Rep 2
st_val = rep2_rng(1);
end_val = rep2_rng(2);

frames_rep1 = f_data(st_val:end_val);

r1_min = islocalmin(frames_rep1);
r1_min_vals = find(r1_min==1);
r1_st = [st_val, r1_min_vals+st_val, end_val];

r1_max = islocalmax(frames_rep1);
r1_max_vals = find(r1_max==1);
r1_nd = r1_max_vals+st_val;

all_idxs = [r1_st, r1_nd];
idxs2 = sort(all_idxs);

%% Rep 1 
st_val = rep3_rng(1);
end_val = rep3_rng(2);

frames_rep1 = f_data(st_val:end_val);

r1_min = islocalmin(frames_rep1);
r1_min_vals = find(r1_min==1);
r1_st = [st_val, r1_min_vals+st_val, end_val];

r1_max = islocalmax(frames_rep1);
r1_max_vals = find(r1_max==1);
r1_nd = r1_max_vals+st_val;

all_idxs = [r1_st, r1_nd];
idxs3 = sort(all_idxs);

%% Combine the data into one data structure. 

for i = 1:48 
    data{i, 1} = v_data(idxs(i):idxs(i+1)-1);
    data{i, 2} = v_data(idxs2(i):idxs2(i+1)-1);
    data{i, 3} = v_data(idxs3(i):idxs3(i+1)-1);
end 

% Add 4th column for mean.

for j = 1:height(data)

    d1 = data{j, 1};
    n_per_col(1) = numel(d1);
    d2 = data{j, 2};
    n_per_col(2) = numel(d2);
    d3 = data{j, 3};
    n_per_col(3) = numel(d3);

    min_val = min(n_per_col);

    d = vertcat(d1(1:min_val), d2(1:min_val), d3(1:min_val));
    mean_resp = nanmean(d);

    data{j, 4} = mean_resp;
end 

%% PLOT

% Number of subplots
numPlots = 16;
theta = linspace(0, 2*pi, numPlots+1); % Angles (add 2*pi to complete the circle)
theta = theta(1:end-1); % Remove the last point since it overlaps with the first

% Center and radius of the circle
centerX = 0.5; % Normalized X-center of the circle
centerY = 0.5; % Normalized Y-center of the circle
radius = 0.35; % Normalized radius of the circle

plot_order= [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

%% Create the figure

max_v = zeros(numPlots, 3);

for sp = 1:3

    figure;

    % Loop to create subplots
    for i = 1:numPlots

        % Compute subplot center position in normalized units
        x = centerX + radius * cos(theta(i)); % X-coordinate
        y = centerY + radius * sin(theta(i)); % Y-coordinate
        
        % Define the size of each subplot
        subplotWidth = 0.15; % Normalized width
        subplotHeight = 0.15; % Normalized height
        
        % Define position of subplot [x, y, width, height]
        subplotPosition = [x - subplotWidth/2, y - subplotHeight/2, subplotWidth, subplotHeight];
        
        % Create the axes at the specified position
        ax = axes('Position', subplotPosition);
    
        % plot data for each rep
        d_idx = plot_order(i) + 16*(sp-1);
    
        for r = 1:4
            d2plot = data{d_idx, r}*10;
            x_vals = 1:numel(d2plot);
    
            if r ==1 
                % Plot median voltage acrss entire recording in background. 
                plot(ax, [1 x_vals(end)], [median_voltage, median_voltage], 'Color', [0.7 0.7 0.7], 'LineWidth', 1) 
                hold on
            end 
    
            if r<4
                plot(ax, x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            else 
                plot(ax, x_vals, d2plot, 'Color', 'k', 'LineWidth', 1.2);
            end
    
        end 
    
        ylim([-80 -10])

        max_rep_voltage = max(data{d_idx, 4}*10);
        max_v(i, sp) = max_rep_voltage;

        % Add a text annotation in the bottom-left corner
        text(ax, x_vals(end)*0.75, median_voltage*1.2, sprintf('%.2f', max_rep_voltage), 'FontSize', 8, 'Color', 'k');
        
        % Turn off the axes for better visualization
        axis(ax, 'off');
    end

    if sp == 1
        sgtitle('14 dps - 4 pixel bar stimuli - 45 pix square')
    elseif sp == 2
        sgtitle('28 dps - 4 pixel bar stimuli - 45 pix square')
    elseif sp == 3 
        sgtitle('56 dps - 4 pixel bar stimuli - 45 pix square')
    end 
    
    f = gcf;
    f.Position = [303 78 961 969];

end 


%% Plot heat map of the max voltage reached during each rep. 

figure; imagesc(max_v); hcb = colorbar;
cm_inferno=inferno(1000);
colormap(cm_inferno)
ax_c= gca;
ax_c.TickDir = 'out';
ax_c.LineWidth = 1;
ax_c.FontSize = 12; 
box off

yticks([1,5,9, 13])
yticklabels({'0', '90', '180', '270'});
ylabel('Direction - deg')

xticks([1,2,3])
xticklabels({'14', '28', '56'})
xlabel('Speed - dps')

colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String','Max voltage (mV)');

f2 = gcf;
f2.Position = [620   386   190   581];






























