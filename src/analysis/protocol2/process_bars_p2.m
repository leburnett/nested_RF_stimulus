function resultant_angle = process_bars_p2(exp_folder, metadata, PROJECT_ROOT)
% PROCESS_BARS_P2  Analyze bar stimulus responses for direction selectivity.
%
%   RESULTANT_ANGLE = PROCESS_BARS_P2(EXP_FOLDER, METADATA, PROJECT_ROOT)
%   extracts responses to moving bar stimuli, calculates direction
%   selectivity metrics, and generates visualization figures.
%
%   INPUTS:
%     exp_folder   - Full path to the Protocol 2 experiment directory
%     metadata     - Structure containing: .Strain, .Frame, .Side, .Age
%     PROJECT_ROOT - Base directory for saving results and figures
%
%   OUTPUT:
%     resultant_angle - Preferred direction in radians (vector sum result)
%
%   ANALYSIS PIPELINE:
%     1. Loads voltage and frame data from TDMS log files
%     2. Parses individual bar stimulus responses using PARSE_BAR_DATA
%     3. Calculates max/min voltage for each of 16 directions
%     4. Computes direction selectivity metrics:
%        - DSI (vector sum and PD/ND methods)
%        - FWHM of tuning curve
%        - Circular variance
%        - Symmetry ratio
%        - Von Mises parameters (theta, kappa)
%     5. Generates visualization figures:
%        - Timeseries with polar plot overlay
%        - Polar plot with vector sum arrow
%        - Heatmap of direction responses
%
%   OUTPUT FILES:
%     Saves to PROJECT_ROOT/results/bar_results/:
%       peak_vals_<strain>_<on_off>_<date>_<time>.mat
%     Contains: bar_results structure, raw data, aligned data, PD ordering
%
%     Saves to PROJECT_ROOT/figures/bar_stimuli/:
%       PDF figures of timeseries, polar plots, and heatmaps
%
%   COMPUTED METRICS:
%     For both slow (28 dps) and fast (56 dps) bars:
%       - max_v_polar: Max voltage per direction
%       - magnitude: Vector sum magnitude
%       - angle_rad: Preferred direction (radians)
%       - fwhm: Full-width half-maximum (degrees)
%       - cv: Circular variance
%       - DSI_vector: Direction selectivity (vector method)
%       - DSI_pdnd: Direction selectivity (PD/ND method)
%       - sym_ratio: Tuning symmetry
%
%   See also PARSE_BAR_DATA, FIND_PD_AND_ORDER_IDX, COMPUTE_BAR_RESPONSE_METRICS

    results_folder = fullfile(PROJECT_ROOT, "results", "bar_results");
    if ~isfolder(results_folder)
        mkdir(results_folder);
    end
    
    figures_folder = fullfile(PROJECT_ROOT, "figures", "bar_stimuli");
    if ~isfolder(figures_folder)
        mkdir(figures_folder);
    end
    
    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
    on_off = params.on_off;
    params.date = date_str;
    params.time = time_str;
    params.strain = metadata.Strain;
    
    f_data = Log.ADC.Volts(1, :); % frame data
    
    v_data = Log.ADC.Volts(2, :)*10; % voltage data
    median_v = median(v_data);

    % Display figure of the frame position and voltage data to get a quick
    % overview of the quality of the recording:
    figure; 
    subplot(2,1,1)
    plot(f_data);
    xlim([0 numel(f_data)])
    ylabel('Frame #')
    subplot(2,1,2)
    plot(v_data);
    xlim([0 numel(f_data)])
    ylabel('Voltage')
    f = gcf;
    f.Position = [50 633 1679 376];

    % v2_data = v_data - median_v; % Get the median-subtracted voltage.
    
    % Parse the raw voltage data for the different bar stimuli
    data = parse_bar_data(f_data, v_data, on_off); % Update this to work for both slow and fast bars
    
    %% PLOTTING THE BAR DATA

    % Plot the timeseries responses with a polar plot in the middle.
    save_fig = 1;
    [max_v, min_v] = plot_timeseries_polar_bars(data, median_v, params, save_fig, figures_folder);
    f = gcf;
    f.Position = [303   380   688   667];

    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars

    % Plot only the polar plot with an arrow overlaid.
    resultant_angle = plot_polar_with_arrow(max_v, median_v, params, save_fig, figures_folder);
    
    % Plot a heat map of the maximum responses to the bars moving in the 16
    % directions at 2 different speeds.
    plot_heatmap_bars(max_v)
    
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
    [d_fast, ord_fast, magnitude_fast, angle_rad_fast, fwhm_fast, cv_fast, thetahat_fast, kappa_fast] = find_PD_and_order_idx(max_v_polar2, median_v);
    
    [sym_ratio_fast, DSI_fast, DSI_pdnd_fast, vector_sum_fast] = compute_bar_response_metrics(d_fast);
    
    % % % % % Save the data:
    bar_results = struct();
    
    bar_results.Date = date_str;
    bar_results.Time = time_str;
    bar_results.Strain = metadata.Strain;
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
    bar_results.resultant_angle = resultant_angle;
    
    save(fullfile(results_folder, strcat('peak_vals_', metadata.Strain, '_', on_off, '_', date_str, '_', time_str, '.mat'))...
        , "bar_results"...
        , 'data' ...
        , 'data_aligned'...
        , 'ord'...
        , 'd_slow'...
        );

end