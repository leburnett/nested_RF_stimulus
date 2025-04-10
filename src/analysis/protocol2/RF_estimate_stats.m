
% Statistics - comparing the size of the RF estimations between groups and
% types of cell. 

data_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/RF_estimate';
cd(data_folder)
files = dir('RF_est*');

var_names_types = [["date", "string"]...
    , ["time", "string"]...
    , ["strain", "double"]...
    , ["type", "double"]...
    , ["peak_v", "double"]...
    , ["x", "double"]...
    , ["y", "double"]...
    , ["sigma_x", "double"]...
    , ["sigma_y", "double"]...
    , ["baseline_v", "double"]...
    , ["gp", "double"]...
    ];

data_tbl = table('Size',[0,size(var_names_types,1)],... 
	'VariableNames', var_names_types(:,1),...
	'VariableTypes', var_names_types(:,2));

for i = 1:height(files)
    fname = files(i).name;
    load(fname, 'optEx')

    str_parts = split(fname(1:end-4), '_');

    data_tbl.date{i} = strcat(str_parts{3:5});
    data_tbl.time{i} = strcat(str_parts{6:7});

    if str_parts{end} == "OFF"
        data_tbl.type(i) = 0; % OFF
    else
        data_tbl.type(i) = 1; % ON
    end 

    if str_parts{end-1} == "control"
        data_tbl.strain(i) = 0; % control
    else
        data_tbl.strain(i) = 1; % ttl
    end 

    if str_parts{end-1} == "control"
        if str_parts{end} == "OFF"
            data_tbl.gp(i) = 1; % control OFF
        else
            data_tbl.gp(i) = 2; % control ON
        end 
    else
        if str_parts{end} == "OFF"
            data_tbl.gp(i) = 3; % ttl OFF
        else
            data_tbl.gp(i) = 4; % ttl ON
        end 
    end 

    data_tbl.peak_v(i) = optEx(1);
    data_tbl.x(i) = optEx(2);
    data_tbl.y(i) = optEx(3);
    data_tbl.sigma_x(i) = optEx(4);
    data_tbl.sigma_y(i) = optEx(5);
    data_tbl.baseline_v(i) = optEx(6);
end 


data_type = 'peak_v';

%% Unbalanced n-way ANOVA:

ttl_on = find(data_tbl.type==1 & data_tbl.strain==1);
ttl_off = find(data_tbl.type==0 & data_tbl.strain==1);

ctl_on = find(data_tbl.type==1 & data_tbl.strain==0);
ctl_off = find(data_tbl.type==0 & data_tbl.strain==0);

mean_on_ttl = mean(data_tbl.(data_type)(ttl_on)); 
mean_off_ttl = mean(data_tbl.(data_type)(ttl_off)); 
mean_on_ctl = mean(data_tbl.(data_type)(ctl_on)); 
mean_off_ctl = mean(data_tbl.(data_type)(ctl_off)); 

[p, tbl, stats] = anovan(data_tbl.(data_type), {data_tbl.type, data_tbl.strain}, ...
                         'model', 'interaction', ...
                         'varnames', {'Type', 'Strain'});

c = multcompare(stats, 'Dimension', [1,2]);

% boxplot(data_tbl.sigma_x, {data_tbl.type, data_tbl.strain

%% Plot figure with boxchart of data

figure
boxchart(data_tbl.gp, data_tbl.(data_type))
xticks([1,2,3,4])
xticklabels({'ctrl-off', 'ctrl-on', 'ttl-off', 'ttl-on'})
xlabel('group')
ylabel(strrep(data_type, '_', '-'))

ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.LineWidth = 1.2;
ax.FontSize = 12;

f = gcf;
f.Position = [798   633   228   401];


