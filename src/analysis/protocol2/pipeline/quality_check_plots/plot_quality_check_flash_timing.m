function plot_quality_check_flash_timing(f_data, flash_dur_ms, save_fig, PROJECT_ROOT)

    figure; 
    plot(f_data);
    hold on;
    plot([flash_dur_ms , flash_dur_ms], [0 400], 'r', 'LineWidth', 1.2)
    
    f=gcf;
    
    if save_fig
        figures_folder = fullfile(PROJECT_ROOT, "figures", "quality");
        fname = fullfile(figures_folder, strcat('Quality_flash_timing_', date_str, '_', time_str, '_', strain_str, '_', type_str,  ".pdf"));
        exportgraphics(f ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    end 

end 
