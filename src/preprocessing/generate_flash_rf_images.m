function success = generate_flash_rf_images(cell_id, exp_folder, rf_results_path, output_path, heatmap_path, resultant_angle)
% GENERATE_FLASH_RF_IMAGES  Generate flash RF timeseries and heatmap PNGs.
%
%   SUCCESS = GENERATE_FLASH_RF_IMAGES(CELL_ID, EXP_FOLDER, RF_RESULTS_PATH,
%       OUTPUT_PATH, HEATMAP_PATH, RESULTANT_ANGLE)
%   loads raw voltage/frame data from the P2 experiment folder and the
%   pre-computed rf_results .mat file, then calls
%   plot_rf_estimate_timeseries_line and plot_heatmap_flash_responses
%   to generate two figures. Saves both as PNG.
%
%   INPUTS:
%     cell_id          - Cell identifier string (YYYY_MM_DD_HH_MM)
%     exp_folder       - Full path to the P2 experiment directory
%     rf_results_path  - Full path to the rf_results .mat file
%     output_path      - Full path for the timeseries output PNG file
%     heatmap_path     - Full path for the heatmap output PNG file
%     resultant_angle  - Preferred direction from bar analysis (radians).
%                        Use 0 if not available.
%
%   OUTPUT:
%     success - true if images were generated, false otherwise
%
%   NOTES:
%     This function replicates the data loading and idx computation from
%     process_flash_p2.m (lines 62-104) to avoid running the full pipeline.
%     It requires load_protocol2_data, plot_rf_estimate_timeseries_line,
%     and plot_heatmap_flash_responses on the MATLAB path.

success = false;

try
    % Load raw data from experiment folder
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));
    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);

    f_data = Log.ADC.Volts(1, :); % frame data
    v_data = Log.ADC.Volts(2, :) * 10; % voltage data (mV)
    median_v = median(v_data);
    v2_data = v_data - median_v; % median-subtracted voltage

    % Load rf_results for data_comb and cmap_id
    S = load(rf_results_path, 'rf_results');
    rf = S.rf_results;
    on_off = rf.Type;

    % Compute flash start indices (same logic as process_flash_p2.m:74-92)
    diff_f_data = diff(f_data);
    n_flashes_4px = 196;
    n_flashes_6px = 100;

    if on_off == "off"
        idx = find(diff_f_data == 1 & f_data(2:end) == 1);
        idx = idx([1,2,5,6,9,10]);
    elseif on_off == "on"
        idx_4 = find(diff_f_data == 1 + n_flashes_4px & f_data(2:end) == 1 + n_flashes_4px);
        idx_4(:, [2,4,6]) = [];
        idx_6 = find(diff_f_data == 1 + n_flashes_6px & f_data(2:end) == 1 + n_flashes_6px);
        idx = sort(horzcat(idx_4, idx_6));
    end

    % Extract slow flash data from rf_results
    % NOTE: process_flash_p2 loops over px_size=[4,6] but stores both under
    % rf_results.slow, so the saved file always has the LAST iteration (6px,
    % 10x10).  Auto-detect px_size from the data dimensions.
    data_comb = rf.slow.data_comb{1};
    cmap_id = rf.slow.cmap_id{1};
    data_comb2 = rescale(data_comb, 0, 1);

    grid_sz = size(data_comb, 1);
    if grid_sz == 14
        px_size = 4;
    elseif grid_sz == 10
        px_size = 6;
    else
        error('Unexpected data_comb size %d — expected 14 (4px) or 10 (6px).', grid_sz);
    end

    % Set up params for the plotting function
    plot_params = struct();
    plot_params.date = date_str;
    plot_params.time = time_str;
    plot_params.strain = rf.Strain;
    plot_params.on_off = on_off;
    plot_params.resultant_angle = resultant_angle;

    % Generate the timeseries figure
    f = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, ...
        "slow", px_size, idx, plot_params);

    % Export timeseries as PNG with white background
    f.Color = 'w';
    exportgraphics(f, output_path, 'Resolution', 300, 'BackgroundColor', 'white');
    close(f);
    fprintf('  Flash RF timeseries saved: %s\n', output_path);

    % Generate and export heatmap
    f2 = plot_heatmap_flash_responses(data_comb2);
    f2.Color = 'w';
    exportgraphics(f2, heatmap_path, 'Resolution', 300, 'BackgroundColor', 'white');
    close(f2);
    fprintf('  Flash RF heatmap saved: %s\n', heatmap_path);

    % Clean up any files the plotting function auto-generated in cwd
    % (it saves PDF and PNG in the current directory)

    success = true;

catch ME
    warning('Failed to generate flash RF image for %s: %s', cell_id, ME.message);
end

end
