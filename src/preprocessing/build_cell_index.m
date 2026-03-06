function build_cell_index(cells_table, output_dir)
% BUILD_CELL_INDEX  Write cell_index.json from experiment metadata and image availability.
%
%   BUILD_CELL_INDEX(CELLS_TABLE, OUTPUT_DIR)
%   scans the output directories for available images and writes a JSON
%   index file mapping each cell to its metadata and image paths.
%
%   INPUTS:
%     cells_table - Table from experiment log (filtered, no test cells)
%     output_dir  - Dashboard images root directory
%
%   OUTPUT:
%     Writes cell_index.json to output_dir

n = height(cells_table);
cells = cell(n, 1);

for i = 1:n
    row = cells_table(i, :);
    cell_id = sprintf('%s_%s', row.date{1}, row.time{1});

    c = struct();
    c.cell_id = cell_id;

    % Parse date for display
    date_parts = strsplit(row.date{1}, '_');
    c.date = sprintf('%s-%s-%s', date_parts{1}, date_parts{2}, date_parts{3});
    time_parts = strsplit(row.time{1}, '_');
    c.time = sprintf('%s:%s', time_parts{1}, time_parts{2});

    c.strain = row.strain{1};

    % Map genotype label
    if isfield(row, 'genotype') && ~isempty(row.genotype{1}) && ~ismissing(row.genotype{1})
        c.genotype = row.genotype{1};
    else
        c.genotype = '';
    end

    c.preferred_contrast = row.preferred_contrast{1};

    if isfield(row, 'project_type') && ~ismissing(row.project_type{1})
        c.project_type = row.project_type{1};
    else
        c.project_type = '';
    end

    c.age = row.age(1);

    % Sweep inclusion
    if isfield(row, 'Incl_SweepAnalysis_') && ~ismissing(row.Incl_SweepAnalysis_{1})
        c.incl_sweep = row.Incl_SweepAnalysis_{1};
    else
        c.incl_sweep = '';
    end

    % Bar flash inclusion
    if isfield(row, 'Incl_BarFlashAnalysis_') && ~ismissing(row.Incl_BarFlashAnalysis_{1})
        c.incl_bar_flash = row.Incl_BarFlashAnalysis_{1};
    else
        c.incl_bar_flash = '';
    end

    % Notes
    if isfield(row, 'NotesStart') && ~ismissing(row.NotesStart{1})
        c.notes_start = row.NotesStart{1};
    else
        c.notes_start = '';
    end
    if isfield(row, 'NotesEnd') && ~ismissing(row.NotesEnd{1})
        c.notes_end = row.NotesEnd{1};
    else
        c.notes_end = '';
    end

    % Check image availability — gridplots (4 per cell)
    for gi = 1:4
        gp_fname = sprintf('%s_gridplot_%d.png', cell_id, gi);
        c.(sprintf('has_gridplot_%d', gi)) = isfile(fullfile(output_dir, 'gridplots', gp_fname));
        c.(sprintf('gridplot_%d_path', gi)) = ['gridplots/' gp_fname];
        c.(sprintf('gridplot_%d_thumb', gi)) = ['thumbs/gridplots/' gp_fname];
    end

    % Check image availability — other plot types
    c.has_bar_polar = isfile(fullfile(output_dir, 'bar_polar', [cell_id '_bar_polar.png']));
    c.has_flash_rf = isfile(fullfile(output_dir, 'flash_rf', [cell_id '_flash_rf.png']));
    c.has_flash_rf_heatmap = isfile(fullfile(output_dir, 'flash_rf_heatmap', [cell_id '_flash_rf_heatmap.png']));
    c.has_bar_flash = isfile(fullfile(output_dir, 'bar_flash', [cell_id '_bar_flash.png']));

    % Image paths (relative)
    c.bar_polar_path = ['bar_polar/' cell_id '_bar_polar.png'];
    c.flash_rf_path = ['flash_rf/' cell_id '_flash_rf.png'];
    c.flash_rf_heatmap_path = ['flash_rf_heatmap/' cell_id '_flash_rf_heatmap.png'];
    c.bar_flash_path = ['bar_flash/' cell_id '_bar_flash.png'];

    % Thumb paths
    c.bar_polar_thumb = ['thumbs/bar_polar/' cell_id '_bar_polar.png'];
    c.flash_rf_thumb = ['thumbs/flash_rf/' cell_id '_flash_rf.png'];
    c.flash_rf_heatmap_thumb = ['thumbs/flash_rf_heatmap/' cell_id '_flash_rf_heatmap.png'];
    c.bar_flash_thumb = ['thumbs/bar_flash/' cell_id '_bar_flash.png'];

    % Display label for the dashboard
    c.display_label = sprintf('%s-%s-%s-%s', c.strain, c.preferred_contrast, row.date{1}, row.time{1});

    cells{i} = c;
end

% Build JSON structure
index = struct();
index.generated_at = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
index.images_dir = output_dir;
index.cells = cells;

json_str = jsonencode(index, 'PrettyPrint', true);
json_path = fullfile(output_dir, 'cell_index.json');
fid = fopen(json_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Cell index written to %s (%d cells)\n', json_path, n);

end
