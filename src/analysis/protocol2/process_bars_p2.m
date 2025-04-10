function process_bars_p2(exp_folder, PROJECT_ROOT)

    results_folder = fullfile(PROJECT_ROOT, "results", "bar_results");
    if ~isfolder(results_folder)
        mkdir(results_folder);
    end
    
    figures_folder = fullfile(PROJECT_ROOT, "figures", "bar_stimuli");
    if ~isfolder(figures_folder)
        mkdir(figures_folder);
    end
    
    strain_str = "T4T5"; % eventually will get this from metadata. 
    
    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
    on_off = params.on_off;
    params.date = date_str;
    params.time = time_str;
    params.strain = strain_str;
    
    f_data = Log.ADC.Volts(1, :); % frame data
    
    v_data = Log.ADC.Volts(2, :)*10; % voltage data
    median_v = median(v_data);
    % v2_data = v_data - median_v; % Get the median-subtracted voltage.
    
    % Parse the raw voltage data for the different bar stimuli
    data = parse_bar_data(f_data, v_data); % Update this to work for both slow and fast bars
    
    % Plot the timeseries responses with a polar plot in the middle.
    save_fig = 0;
    [max_v, min_v] = plot_timeseries_polar_bars(data, median_v, params, save_fig);
    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars

    % Plot only the polar plot with an arrow overlaid.
    plot_polar_with_arrow(max_v, median_v, params, save_fig)
    
    % Plot a heat map of the maximum responses to the bars moving in the 16
    % directions at 2 different speeds.
    plot_heatmap_bars(max_v)
    
    % Align the data angles to plot a linear timeseries plot:
    data_ordered = align_data_by_seq_angles(data);
    
    % Find PR and the order to re-order the data to align PD to pi/2.
    % 1 - slower bar stimuli
    [d_slow, ord, magnitude_slow, angle_rad_slow, fwhm_slow, cv_slow, thetahat_slow, kappa_slow] = find_PD_and_order_idx(max_v_polar1, median_v); % use the max v polar for the slower bars.
    
    [sym_ratio_slow, DSI_slow, DSI_pdnd_slow, vector_sum_slow] = compute_bar_response_metrics(d_slow);
    
    data_aligned = cell(size(data_ordered));
    for kk = 1:32 
        if kk<17
            ord_id = ord(kk);
        else 
            ord_id = ord(kk-16)+16;
        end 
    
        data_aligned(ord_id, :) = data_ordered(kk, :); % Combine the mean timeseries only.
    end 
    
    % 2 - faster bar stimuli
    [d_fast, ord_fast, magnitude_fast, angle_rad_fast, fwhm_fast, cv_fast, thetahat_fast, kappa_fast] = find_PD_and_order_idx(max_v_polar2, median_v);
    
    [sym_ratio_fast, DSI_fast, DSI_pdnd_fast, vector_sum_fast] = compute_bar_response_metrics(d_fast);
    
    % % % % % Save the data:
    bar_results = struct();
    
    bar_results.Date = date_str;
    bar_results.Time = time_str;
    bar_results.Strain = strain_str;
    bar_results.Type = on_off;
    bar_results.slow.max_v_polar = {max_v_polar1};
    bar_results.fast.max_v_polar = {max_v_polar2};
    bar_results.min_v = {min_v};
    bar_results.max_v = {max_v};
    bar_results.median_voltage = median_v;
    
    % output of vector sum:
    bar_results.slow.magnitude = magnitude_slow;
    bar_results.slow.angle_rad = angle_rad_slow;
    bar_results.slow.fwhm = fwhm_slow;
    bar_results.slow.cv = cv_slow;
    bar_results.slow.thetahat = thetahat_slow;
    bar_results.slow.kappa = kappa_slow;
    
    bar_results.fast.magnitude = magnitude_fast;
    bar_results.fast.angle_rad = angle_rad_fast;
    bar_results.fast.fwhm = fwhm_fast;
    bar_results.fast.cv = cv_fast;
    bar_results.fast.thetahat = thetahat_fast;
    bar_results.fast.kappa = kappa_fast;
    
    % symmetry index
    bar_results.slow.sym_ratio = sym_ratio_slow;
    bar_results.fast.sym_ratio = sym_ratio_fast;
    
    % DSI - vector sum 
    bar_results.slow.vector_sum = vector_sum_slow;
    bar_results.slow.DSI_vector = DSI_slow;
    bar_results.slow.DSI_pdnd = DSI_pdnd_slow;
    bar_results.fast.vector_sum = vector_sum_fast;
    bar_results.fast.DSI_vector = DSI_fast;
    bar_results.fast.DSI_pdnd = DSI_pdnd_fast;
    
    save(fullfile(results_folder, strcat('peak_vals_', strain_str, '_', on_off, '_', date_str, '_', time_str, '.mat'))...
        , "bar_results"...
        , 'data' ...
        , 'data_ordered'...
        , 'data_aligned'...
        , 'ord'...
        , 'd_slow'...
        );

end