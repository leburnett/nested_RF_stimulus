
close all
clear

G4_file = dir('**/G4*');
load(fullfile(G4_file.folder, G4_file.name))

% Voltage
v_data = Log.ADC.Volts(2, :)*10;

figure; plot(v_data)
xlim([0 numel(v_data)])
box off
f=gcf;
f.Position = [86  283  1675  166];

figure; histogram(v_data,'Normalization','pdf', 'BinEdges', -75:0.1:-40, 'FaceColor', [0.8 0.8 0.8])
ylabel('PDF')
xlabel('Vm (mV)')

hold on;

% Compute statistics
mean_vm = mean(v_data);
disp(mean_vm)
std_vdata = std(v_data);
disp(std_vdata)

% Define positions for the lines
x_mean = mean_vm;
x_plus = mean_vm + std_vdata;
x_minus = mean_vm - std_vdata;

% Add vertical lines
yl = ylim; % get y-axis limits for line placement
plot([x_mean x_mean], yl, 'r-', 'LineWidth', 1.5);
plot([x_plus x_plus], yl, 'm--', 'LineWidth', 1.2);
plot([x_minus x_minus], yl, 'm--', 'LineWidth', 1.2);

% Add text labels
text(x_mean, yl(2)*0.99, sprintf('Mean = %.2f', x_mean), 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight','bold');
text(x_mean, yl(2)*0.96, sprintf('SD = %.2f', std_vdata), 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight','bold');
text(x_plus+2, yl(2)*0.8, sprintf('%.2f', x_plus), 'Color', 'm', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(x_minus-2, yl(2)*0.8, sprintf('%.2f', x_minus), 'Color', 'm', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

box off
hold off;

