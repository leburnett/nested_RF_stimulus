function plot_quality_check_data_bounds(v_data, median_v, upper_bound_med, lower_bound_med)

figure; 
plot(v_data);
hold on;
plot([1, numel(v_data)],[median_v, median_v], 'r'); % median
plot([1, numel(v_data)],[median_v+upper_bound_med, median_v+upper_bound_med], 'm'); % upper
plot([1, numel(v_data)],[median_v+lower_bound_med, median_v+lower_bound_med], 'm'); % lower

end 