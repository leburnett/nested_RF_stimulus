function plot_quality_check_data_full_rec(f_data, v_data, save_fig, PROJECT_ROOT)

figure; 

subplot(5,1,1)
plot(f_data)
ax = gca;
ax.XAxis.Visible = 'off';
title('f-data')
xlim([0 numel(f_data)])

subplot(5,1,2:3)
plot(v_data)
xlim([0 numel(v_data)])
ylabel('v-data')
ax = gca;
ax.XAxis.Visible = 'off';
title(strcat("var - v-data = ", string(var(v_data))))

subplot(5,1,4:5)
plot(movmean(v_data, 20000))
xlim([0 numel(v_data)])
ylabel('movmean(v-data, 20000)')
hold on 
plot(movmean(v_data, 100000))

title(strcat("var - movmean - 100000 = ", string(var(movmean(v_data, 100000)))))
f = gcf;
f.Position = [25 629 1121 392];

if save_fig
    figures_folder = fullfile(PROJECT_ROOT, "figures", "quality");
    fname = fullfile(figures_folder, strcat('Quality_full_rec_', date_str, '_', time_str, '_', strain_str, '_', type_str,  ".pdf"));
    exportgraphics(f ...
        , fname ...
        , 'ContentType', 'vector' ...
        , 'BackgroundColor', 'none' ...
        ); 
end 

end 













