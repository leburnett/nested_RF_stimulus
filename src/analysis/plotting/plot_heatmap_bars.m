function plot_heatmap_bars(max_v)
% PLOT_HEATMAP_BARS  Create heatmap of direction tuning at two speeds.
%
%   PLOT_HEATMAP_BARS(MAX_V) generates a heatmap visualization showing
%   the maximum neural response for each of 16 bar directions at two
%   different speeds.
%
%   INPUT:
%     max_v - 16x2 array of maximum voltage responses
%             Column 1: slow (28 dps) bar responses
%             Column 2: fast (56 dps) bar responses
%             Rows correspond to 0, 22.5, 45, ... 337.5 degrees
%
%   FIGURE:
%     - Y-axis: Direction in degrees (0, 90, 180, 270 labeled)
%     - X-axis: Speed in degrees per second (28, 56)
%     - Color: Maximum voltage (mV), using inferno colormap
%
%   PURPOSE:
%     Provides a quick visual comparison of direction tuning at
%     different speeds. Can reveal speed-dependent changes in
%     directional preference or selectivity strength.
%
%   See also PLOT_TIMESERIES_POLAR_BARS, PROCESS_BARS_P2, INFERNO

figure; 
imagesc(max_v); 

hcb = colorbar;
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

xticks([1,2, 3])
xticklabels({'28', '56', '168'})
xlabel('Speed - dps')

colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String','Max voltage (mV)');

f2 = gcf;
f2.Position = [620   386   190   581];

end 