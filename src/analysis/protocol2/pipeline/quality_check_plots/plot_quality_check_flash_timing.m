function plot_quality_check_flash_timing(f_data, flash_dur_ms)

figure; 
plot(f_data);
hold on;
plot([flash_dur_ms , flash_dur_ms], [0 400], 'r', 'LineWidth', 1.2)

end 
