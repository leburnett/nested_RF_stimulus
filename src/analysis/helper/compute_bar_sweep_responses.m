function max_v = compute_bar_sweep_responses(bar_data, plot_order, opts)
% COMPUTE_BAR_SWEEP_RESPONSES  Extract peak depolarization from bar sweep means.
%
%   MAX_V = COMPUTE_BAR_SWEEP_RESPONSES(BAR_DATA, PLOT_ORDER, OPTS)
%   computes the depolarization amplitude for each of 16 directions by
%   comparing the 98th percentile voltage during the stimulus period to
%   the mean baseline voltage before stimulus onset.
%
%   INPUTS:
%     bar_data   - Nx(R+1) cell array from parse_bar_data. Rows 1-16 are
%                  slow bar directions. Columns 1-R are individual reps,
%                  column R+1 is the mean across reps.
%     plot_order - 1x16 index mapping from subplot position to data row
%     opts       - Structure with optional fields:
%                    .baseline_range - [start end] sample indices for the
%                                      pre-stimulus baseline period
%                                      (default: [1000 9000])
%                    .stim_trim_end  - Number of samples to trim from the
%                                      end of the stimulus period
%                                      (default: 7000)
%                    .percentile     - Percentile for peak voltage detection
%                                      (default: 98)
%
%   OUTPUT:
%     max_v - 16x1 array of depolarization amplitudes, ordered by subplot
%             position (i.e., max_v(i) corresponds to data row plot_order(i))
%
%   COMPUTATION:
%     For each direction's mean voltage trace:
%       baseline = mean(trace(baseline_range(1):baseline_range(2)))
%       peak     = prctile(trace(9000:end-stim_trim_end), percentile)
%       max_v(i) = abs(peak - baseline)
%
%   See also PARSE_BAR_DATA, PLOT_SLOW_BAR_SWEEP_POLAR, FIND_PD_FROM_LUT

    % Set defaults
    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'baseline_range'), opts.baseline_range = [1000 9000]; end
    if ~isfield(opts, 'stim_trim_end'),  opts.stim_trim_end  = 7000; end
    if ~isfield(opts, 'percentile'),     opts.percentile     = 98; end

    n_dir = numel(plot_order);
    n_reps = size(bar_data, 2) - 1;
    max_v = zeros(n_dir, 1);

    bl_start = opts.baseline_range(1);
    bl_end   = opts.baseline_range(2);

    for subplot_idx = 1:n_dir
        data_row = plot_order(subplot_idx);
        d_mean = bar_data{data_row, n_reps + 1};

        % Baseline: mean voltage before stimulus onset
        mean_before = mean(d_mean(bl_start:bl_end));

        % Stimulus period: from baseline end to trimmed end
        d_stim = d_mean(bl_end:end - opts.stim_trim_end);

        % Depolarization = peak during stimulus minus baseline
        max_v(subplot_idx) = abs(diff([prctile(d_stim, opts.percentile), mean_before]));
    end

end
