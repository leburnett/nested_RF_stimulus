function quality = quality_check_recording(exp_folder, opts)
% QUALITY_CHECK_RECORDING  Assess recording quality for one experiment.
%
%   QUALITY = QUALITY_CHECK_RECORDING(EXP_FOLDER) loads the G4_TDMS_Logs
%   .mat file from EXP_FOLDER, computes recording quality metrics,
%   generates a quality check figure, and returns a struct of metrics
%   including a pass/fail FLAG.
%
%   QUALITY = QUALITY_CHECK_RECORDING(EXP_FOLDER, OPTS) uses the options
%   structure to override default thresholds.
%
%   INPUTS:
%     exp_folder - Full path to the experiment directory containing
%                  'Log Files/G4_TDMS_Log*.mat' and 'currentExp.mat'
%     opts       - (Optional) structure with fields:
%                    .save_fig    - Save quality figure (default: true)
%                    .fig_folder  - Where to save figure (default: exp_folder)
%                    .update_mat  - Update currentExp.mat (default: true)
%                    .mean_range  - [low high] acceptable mean mV
%                                   (default: [-70 -35])
%                    .max_std     - Max acceptable std in mV (default: 20)
%                    .min_std     - Min acceptable std in mV (default: 0.3)
%                    .max_drift   - Max abs drift in mV (default: 8)
%                    .max_range   - Max p98-p2 range in mV (default: 50)
%
%   OUTPUT:
%     quality - Structure with fields:
%       .FLAG          - 1 if recording FAILS quality check, 0 if passes
%       .mean_v        - Mean voltage (mV)
%       .std_v         - Standard deviation (mV)
%       .median_v      - Median voltage (mV)
%       .mad_v         - Median absolute deviation (mV)
%       .p2_v          - 2nd percentile voltage
%       .p98_v         - 98th percentile voltage
%       .range_v       - p98 - p2 (mV)
%       .drift_v       - Mean of last 10% minus mean of first 10% (mV)
%       .var_v         - Variance of voltage
%       .fail_reasons  - Cell array of strings listing which checks failed
%
%   See also RUN_QUALITY_CHECK_1DRF, PROCESS_PROTOCOL2, LOAD_PROTOCOL2_DATA

    if nargin < 2, opts = struct(); end
    opts = set_quality_defaults(opts, exp_folder);

    %% Load voltage data
    log_folder = fullfile(exp_folder, 'Log Files');
    log_file = dir(fullfile(log_folder, 'G4_TDMS*'));
    if isempty(log_file)
        error('quality_check_recording:noLogFile', ...
            'No G4_TDMS log file found in %s', log_folder);
    end
    S = load(fullfile(log_folder, log_file(1).name), 'Log');
    v_data = S.Log.ADC.Volts(2, :) * 10;

    %% Extract date and time from folder name
    if contains(exp_folder, '/')
        parts = split(exp_folder, '/');
    elseif contains(exp_folder, '\')
        parts = split(exp_folder, '\');
    end
    exp_name = parts{end};
    date_str = exp_name(1:10);
    time_str = exp_name(12:16);

    %% Compute quality metrics
    n_samples = numel(v_data);
    pct_10 = round(0.10 * n_samples);

    quality.mean_v   = mean(v_data);
    quality.std_v    = std(v_data);
    quality.median_v = median(v_data);
    quality.mad_v    = mad(v_data, 1);
    quality.var_v    = var(v_data);
    quality.p2_v     = prctile(v_data, 2);
    quality.p98_v    = prctile(v_data, 98);
    quality.range_v  = quality.p98_v - quality.p2_v;
    quality.drift_v  = mean(v_data(end-pct_10+1:end)) - mean(v_data(1:pct_10));

    %% Apply quality checks
    fail_reasons = {};

    if quality.mean_v < opts.mean_range(1) || quality.mean_v > opts.mean_range(2)
        fail_reasons{end+1} = sprintf('mean_v=%.1f outside [%.0f, %.0f]', ...
            quality.mean_v, opts.mean_range(1), opts.mean_range(2));
    end

    if quality.std_v > opts.max_std
        fail_reasons{end+1} = sprintf('std_v=%.1f > %.1f (too noisy)', ...
            quality.std_v, opts.max_std);
    end

    if quality.std_v < opts.min_std
        fail_reasons{end+1} = sprintf('std_v=%.2f < %.1f (no activity)', ...
            quality.std_v, opts.min_std);
    end

    if abs(quality.drift_v) > opts.max_drift
        fail_reasons{end+1} = sprintf('|drift|=%.1f > %.1f', ...
            abs(quality.drift_v), opts.max_drift);
    end

    if quality.range_v > opts.max_range
        fail_reasons{end+1} = sprintf('range=%.1f > %.1f (artifacts)', ...
            quality.range_v, opts.max_range);
    end

    quality.fail_reasons = fail_reasons;
    quality.FLAG = double(~isempty(fail_reasons));

    %% Generate quality check figure
    fig = figure('Position', [50 400 1800 250], 'Visible', 'off');

    % Left panel: voltage time series
    ax1 = subplot(1, 5, 1:4);
    t_sec = (0:n_samples-1) / 10000;
    plot(ax1, t_sec, v_data, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.3);
    hold(ax1, 'on');
    plot(ax1, [t_sec(1) t_sec(end)], [quality.mean_v quality.mean_v], ...
        'r-', 'LineWidth', 1.2);
    hold(ax1, 'off');
    xlim(ax1, [0 t_sec(end)]);
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Vm (mV)');
    box(ax1, 'off');

    % Right panel: voltage histogram
    ax2 = subplot(1, 5, 5);
    bin_lo = floor(min(v_data));
    bin_hi = ceil(max(v_data));
    histogram(ax2, v_data, 'Normalization', 'pdf', ...
        'BinEdges', bin_lo:0.2:bin_hi, ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    hold(ax2, 'on');
    yl = ylim(ax2);
    plot(ax2, [quality.mean_v quality.mean_v], yl, 'r-', 'LineWidth', 1.5);
    plot(ax2, [quality.mean_v + quality.std_v, quality.mean_v + quality.std_v], ...
        yl, 'm--', 'LineWidth', 1);
    plot(ax2, [quality.mean_v - quality.std_v, quality.mean_v - quality.std_v], ...
        yl, 'm--', 'LineWidth', 1);
    hold(ax2, 'off');
    xlabel(ax2, 'Vm (mV)');
    ylabel(ax2, 'PDF');
    box(ax2, 'off');

    % Title with metrics and pass/fail status
    if quality.FLAG
        status_str = 'FAIL';
        title_color = [0.8 0 0];
    else
        status_str = 'PASS';
        title_color = [0 0.5 0];
    end

    sgtitle(sprintf('%s %s | mean=%.1f mV | std=%.1f mV | drift=%.1f mV | %s', ...
        strrep(date_str, '_', '-'), strrep(time_str, '_', ':'), ...
        quality.mean_v, quality.std_v, quality.drift_v, status_str), ...
        'Color', title_color, 'FontWeight', 'bold');

    %% Save figure
    if opts.save_fig
        fig_filename = fullfile(opts.fig_folder, ...
            sprintf('quality_check_%s_%s.png', date_str, time_str));
        exportgraphics(fig, fig_filename, ...
            'ContentType', 'image', 'Resolution', 300);
    end
    close(fig);

    %% Update currentExp.mat with quality metadata
    if opts.update_mat
        mat_path = fullfile(exp_folder, 'currentExp.mat');
        if isfile(mat_path)
            S_exp = load(mat_path, 'metadata');
            S_exp.metadata.quality_FLAG         = quality.FLAG;
            S_exp.metadata.quality_mean_v       = quality.mean_v;
            S_exp.metadata.quality_std_v        = quality.std_v;
            S_exp.metadata.quality_median_v     = quality.median_v;
            S_exp.metadata.quality_mad_v        = quality.mad_v;
            S_exp.metadata.quality_drift_v      = quality.drift_v;
            S_exp.metadata.quality_range_v      = quality.range_v;
            S_exp.metadata.quality_p2_v         = quality.p2_v;
            S_exp.metadata.quality_p98_v        = quality.p98_v;
            S_exp.metadata.quality_fail_reasons = quality.fail_reasons;
            metadata = S_exp.metadata; %#ok<NASGU>
            save(mat_path, 'metadata', '-append');
        end
    end

    fprintf('  Quality check: %s (mean=%.1f, std=%.1f, drift=%.1f)\n', ...
        status_str, quality.mean_v, quality.std_v, quality.drift_v);

end


%% ========================= Default Options ==============================

function opts = set_quality_defaults(opts, exp_folder)
% SET_QUALITY_DEFAULTS  Fill in default values for quality check options.

    if ~isfield(opts, 'save_fig'),    opts.save_fig = true;         end
    if ~isfield(opts, 'fig_folder'),  opts.fig_folder = exp_folder; end
    if ~isfield(opts, 'update_mat'),  opts.update_mat = true;       end
    if ~isfield(opts, 'mean_range'),  opts.mean_range = [-70 -35];  end
    if ~isfield(opts, 'max_std'),     opts.max_std = 20;            end
    if ~isfield(opts, 'min_std'),     opts.min_std = 0.3;           end
    if ~isfield(opts, 'max_drift'),   opts.max_drift = 8;           end
    if ~isfield(opts, 'max_range'),   opts.max_range = 50;          end

end
