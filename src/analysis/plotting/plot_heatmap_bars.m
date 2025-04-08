function plot_heatmap_bars(max_v)
% Create heatmap of the maximum response in the 16 directions for the two
% speeds that the bar stimuli are presented.
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

xticks([1,2])
xticklabels({'28', '56'})
xlabel('Speed - dps')

colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String','Max voltage (mV)');

f2 = gcf;
f2.Position = [620   386   190   581];

end 