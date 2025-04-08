function plot_quality_check_var_reps_prctile(var_across_reps, var_within_reps, diff_mean, max_data, min_data, save_fig, PROJECT_ROOT)

    % With median value across all flashes in title.
    figure; 
    subplot(2,3,1) % 1 - variance across reps - consistency
    imagesc(var_across_reps); colorbar
    med_var_X_reps = median(reshape(var_across_reps, [1, 196]));
    title(strcat("CoV across reps - ", string(med_var_X_reps)))
    
    subplot(2,3,4)
    histogram(var_across_reps);
    xlabel('Coeff. of var. across reps - 98% val')
    
    subplot(2,3,2) % 2 - variance within reps - strength of response
    imagesc(var_within_reps); colorbar
    med_var_W_reps = var(reshape(var_within_reps, [1, 196]));
    title(strcat("var within reps - ", string(med_var_W_reps)))
    
    subplot(2,3,5)
    histogram(var_within_reps);
    xlabel('Var. within reps - mean')
    
    subplot(2,3,3) % 3 - difference between max and min per flash for mean.
    imagesc(diff_mean); colorbar
    med_diff_mean = median(reshape(diff_mean, [1, 196]));
    title(strcat("diff mean - ", string(med_diff_mean)))
    
    subplot(2,3,6)
    histogram(diff_mean);
    xlabel('Diff between max and min of mean')
     
    f = gcf;
    f.Position = [1   684   827   363];

    if save_fig
        figures_folder = fullfile(PROJECT_ROOT, "figures", "quality");
        fname = fullfile(figures_folder, strcat('Quality_var_', date_str, '_', time_str, '_', strain_str, '_', type_str,  ".pdf"));
        exportgraphics(f ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    end 

    % % % 

    figure; 
    subplot(2,1,1)
    imagesc(max_data); colorbar; title('max - 98th prctile')
    subplot(2,1,2)
    histogram(max_data)
    f=gcf;
    f.Position = [620   501   275   466];

    if save_fig
        figures_folder = fullfile(PROJECT_ROOT, "figures", "quality");
        fname = fullfile(figures_folder, strcat('Quality_max_heatmap', date_str, '_', time_str, '_', strain_str, '_', type_str,  ".pdf"));
        exportgraphics(f ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    end 

    % % %
    
    figure; 
    subplot(2,1,1)
    imagesc(min_data*-1); colorbar; title('min - 98th prctile')
    subplot(2,1,2)
    histogram(min_data)
    f=gcf;
    f.Position = [620   501   275   466];

    if save_fig
        figures_folder = fullfile(PROJECT_ROOT, "figures", "quality");
        fname = fullfile(figures_folder, strcat('Quality_min_heatmap', date_str, '_', time_str, '_', strain_str, '_', type_str,  ".pdf"));
        exportgraphics(f ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    end 

end 