function selections = get_input_parameters()
% GET_INPUT_PARAMETERS  Display GUI dialog to collect Protocol 2 metadata.
%
%   SELECTIONS = GET_INPUT_PARAMETERS() opens a graphical user interface
%   that prompts the experimenter to enter metadata required for Protocol 2
%   generation. The function blocks (via uiwait) until the user clicks
%   "Confirm".
%
%   OUTPUT:
%     selections - struct with fields:
%       .Frame  - double, peak frame number from Protocol 1. This
%                 identifies the flash location with the strongest response
%                 and is used to centre Protocol 2 stimuli.
%       .Side   - char, arena hemisphere ('L' or 'R').
%                 'L' = left hemisphere, 'R' = right hemisphere.
%       .Age    - char, fly age in days post-eclosion ('1'-'8' or 'NaN').
%       .Strain - char, fly genotype. One of:
%                 '42F06_T4T5_control' - T4/T5 GAL4 control
%                 '42F06_T4T5_ttl'     - T4/T5 GAL4 experimental (TTL)
%                 'ss324_t4t5'         - T4/T5 split-GAL4
%                 'jfrc100_es'         - Empty split control
%                 'ss00395_TmY3'       - TmY3 split-GAL4
%                 'ss03722_Tm5Y'       - Tm5Y split-GAL4
%                 'ss02701_TmY18'      - TmY18 split-GAL4
%                 'ss02594_TmY5a'      - TmY5a split-GAL4
%                 'test'               - For testing purposes
%
%   USAGE:
%     Called automatically by GENERATE_PROTOCOL2. The experimenter must:
%       1. Enter the peak_frame number identified from Protocol 1 responses
%       2. Select the arena side matching the recording hemisphere
%       3. Enter the fly's age
%       4. Select the correct genotype
%       5. Click "Confirm" to proceed
%
%   See also GENERATE_PROTOCOL2, PATT_FRAME_TO_COORD

    % Create a figure window for the input dialog
    fig = uifigure('Name', 'Protocol 2 meta data', 'Position', [500 500 300 300]);

    % Dropdown for Frame number to centre p2 on:
    uilabel(fig, 'Position', [20 250 80 22], 'Text', 'Frame #:');
    frameEditField = uieditfield(fig, 'numeric', ...
    'Position', [120 250 150 22], ...
    'Value', 1);

    % Dropdown for Arena side
    uilabel(fig, 'Position', [20 200 80 22], 'Text', 'Arena side:');
    arenaDropdown = uidropdown(fig, ...
        'Position', [120 200 150 22], ...
        'Items', {'L', 'R'}, ...
        'Value', 'L');

    % Dropdown for Age
    uilabel(fig, 'Position', [20 150 80 22], 'Text', 'Age of Fly:');
    ageDropdown = uidropdown(fig, ...
        'Position', [120 150 150 22], ...
        'Items', {'1', '2', '3', '4', '5', '6', '7', '8', 'NaN'}, ...
        'Value', '4');

    % Dropdown for Strain
    uilabel(fig, 'Position', [20 100 80 22], 'Text', 'Strain:');
    strainDropdown = uidropdown(fig, ...
        'Position', [120 100 150 22], ...
        'Items', {'42F06_T4T5_control', '42F06_T4T5_ttl', 'ss324_t4t5', 'jfrc100_es', 'ss00395_TmY3', 'ss03722_Tm5Y', 'ss02701_TmY18', 'ss02594_TmY5a', 'test'}, ...
        'Value', '42F06_T4T5_control');

    % Variable to hold selections
    selections = struct('Frame', '', 'Side', '', 'Age', '', 'Strain', '');

    % Confirm Button with callback function to retrieve values and close the figure
    confirmButton = uibutton(fig, 'push', ...
        'Position', [100 10 100 30], ...
        'Text', 'Confirm', ...
        'ButtonPushedFcn', @(btn, event) confirmSelections(frameEditField, arenaDropdown, ageDropdown, strainDropdown, fig));

    % Wait for the figure to close before proceeding
    uiwait(fig);

    % Callback function to store selections and close the figure
    function confirmSelections(frameEditField, arenaDropdown, ageDropdown, strainDropdown,fig)
        selections.Frame = frameEditField.Value;
        selections.Side = arenaDropdown.Value;
        selections.Age = ageDropdown.Value;
        selections.Strain = strainDropdown.Value;
        close(fig);
    end

end 









