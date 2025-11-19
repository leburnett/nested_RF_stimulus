function selections = get_input_parameters()
% Input the parameters for the generation of protocol 2:

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
        'Items', {'SS02344_T4','42F06_T4T5_control', '42F06_T4T5_ttl', 'ss324_t4t5', 'jfrc100_es', 'ss00395_TmY3', 'ss03722_Tm5Y', 'ss02701_TmY18', 'ss02594_TmY5a', 'test'}, ...
        'Value', 'SS02344_T4');

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









