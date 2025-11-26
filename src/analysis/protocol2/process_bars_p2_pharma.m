function resultant_angle = process_bars_p2_pharma(exp_folder, metadata, PROJECT_ROOT)

    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
    on_off = params.on_off;
    params.date = date_str;
    params.time = time_str;
    params.strain = metadata.Strain;
    params.drug = metadata.Drug;
    if params.drug == "none"
        params.application = "pre";
        params.drug = "no-drug";
    else
        params.application = "post";
    end
    
    f_data = Log.ADC.Volts(1, :); % frame data
    
    v_data = Log.ADC.Volts(2, :)*10; % voltage data
    median_v = median(v_data);

    results_folder = fullfile(PROJECT_ROOT, "results", "bar_results", params.drug, params.application);
    if ~isfolder(results_folder)
        mkdir(results_folder);
    end
    
    figures_folder = fullfile(PROJECT_ROOT, "figures", "bar_stimuli", params.drug, params.application);
    if ~isfolder(figures_folder)
        mkdir(figures_folder);
    end

    % Display figure of the frame position and voltage data to get a quick
    % overview of the quality of the recording:

    % % Plot figure of voltage and frame data - linked axes.
    % figure; 
    % tiledlayout(2,1);
    % ax1 = nexttile;
    % plot(f_data);
    % xlim([0 numel(f_data)])
    % ylabel('Frame #')
    % ax2 = nexttile;
    % plot(v_data);
    % hold on
    % plot([0 numel(f_data)], [median_v, median_v], 'r')
    % xlim([0 numel(f_data)])
    % ylabel('Voltage')
    % linkaxes([ax1, ax2], 'x');
    % f = gcf;
    % f.Position = [15 621 1769 322];

    % % Plot voltage and frame data in the same figure.
    % figure; 
    % yyaxis left
    % plot(f_data);
    % xlim([0 numel(f_data)])
    % yyaxis right
    % plot(v_data);
    % hold on
    % plot([0 numel(f_data)], [median_v, median_v], 'r')
    % f = gcf;
    % f.Position = [15 621 1769 322];

    % v2_data = v_data - median_v; % Get the median-subtracted voltage.
    
    % Parse the raw voltage data for the different bar stimuli
    data = parse_bar_data_pharma(f_data, v_data);
    
    %% PLOTTING THE BAR DATA

    % Plot the timeseries responses with a polar plot in the middle.
    save_fig = 1;
    [max_v, min_v] = plot_timeseries_polar_bars_pharma(data, median_v, params, save_fig, figures_folder);

    % Plot repetition responses to bar sweeps in grid formation: 
    plot_timeseries_grid_bars_pharma(data, median_v, params, save_fig, figures_folder)

    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % 28 dps
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % 56 dps
    max_v_polar3 = vertcat(max_v(:, 3), max_v(1, 3)); % 168 dps
    max_v_polar4 = vertcat(max_v(:, 4), max_v(1, 4)); % 250 dps
    max_v_polar5 = vertcat(max_v(:, 5), max_v(1, 5)); % 500 dps

    % Plot only the polar plot with an arrow overlaid.
    resultant_angle = plot_polar_with_arrow_pharma(max_v, median_v, params, save_fig, figures_folder);

    % Plot the polar plot of each rep:
    plot_each_rep_28dps_pharma(data, params, save_fig, figures_folder)

    
    %% ANALYSING THE BAR DATA

    % Align the data angles to plot a linear timeseries plot:
    data_ordered = align_data_by_seq_angles(data);
    
    % Find PD and the order to re-order the data to align PD to pi/2.
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
    [d_fast, ~, magnitude_fast, angle_rad_fast, fwhm_fast, cv_fast, thetahat_fast, kappa_fast] = find_PD_and_order_idx(max_v_polar2, median_v);
    
    [sym_ratio_fast, DSI_fast, DSI_pdnd_fast, vector_sum_fast] = compute_bar_response_metrics(d_fast);

    % 3 - even faster bar stimuli
    [d_vfast, ~, magnitude_vfast, angle_rad_vfast, fwhm_vfast, cv_vfast, thetahat_vfast, kappa_vfast] = find_PD_and_order_idx(max_v_polar3, median_v);
    
    [sym_ratio_vfast, DSI_vfast, DSI_pdnd_vfast, vector_sum_vfast] = compute_bar_response_metrics(d_vfast);
    
    % % % % % Save the data:
    bar_results = struct();
    
    bar_results.Date = date_str;
    bar_results.Time = time_str;
    bar_results.Strain = metadata.Strain;
    bar_results.Type = on_off;
    bar_results.slow.max_v_polar = {max_v_polar1};
    bar_results.fast.max_v_polar = {max_v_polar2};
    bar_results.dps168.max_v_polar = {max_v_polar3};
    bar_results.dps250.max_v_polar = {max_v_polar4};
    bar_results.dps500.max_v_polar = {max_v_polar5};
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

    bar_results.dps168.magnitude = magnitude_vfast;
    bar_results.dps168.angle_rad = angle_rad_vfast;
    bar_results.dps168.fwhm = fwhm_vfast;
    bar_results.dps168.cv = cv_vfast;
    bar_results.dps168.thetahat = thetahat_vfast;
    bar_results.dps168.kappa = kappa_vfast;
    
    % symmetry index
    bar_results.slow.sym_ratio = sym_ratio_slow;
    bar_results.fast.sym_ratio = sym_ratio_fast;
    bar_results.dps168.sym_ratio = sym_ratio_vfast;
    
    % DSI - vector sum 
    bar_results.slow.vector_sum = vector_sum_slow;
    bar_results.slow.DSI_vector = DSI_slow;
    bar_results.slow.DSI_pdnd = DSI_pdnd_slow;
    bar_results.fast.vector_sum = vector_sum_fast;
    bar_results.fast.DSI_vector = DSI_fast;
    bar_results.fast.DSI_pdnd = DSI_pdnd_fast;
    bar_results.dps168.vector_sum = vector_sum_vfast;
    bar_results.dps168.DSI_vector = DSI_vfast;
    bar_results.dps168.DSI_pdnd = DSI_pdnd_vfast;
    bar_results.resultant_angle = resultant_angle;
    
    save(fullfile(results_folder, strcat('peak_vals_', metadata.Strain, '_', on_off, '_', date_str, '_', time_str, '.mat'))...
        , "bar_results"...
        , 'data' ...
        , 'data_aligned'...
        , 'ord'...
        , 'd_slow'...
        , 'd_fast'...
        , 'd_vfast'...
        );

end