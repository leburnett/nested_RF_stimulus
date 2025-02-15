%% Analyse the differences in DSI and RF size between the control and turtle cells.

% The results tables are found here:
% '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results';

%% 1 - Differences in DSI to bar stimuli:

% Load the bar results table:
data_tbl = bar_res_table;
% Rename the column variables
data_tbl = renamevars(data_tbl,["Date","Time","Strain","Type"],["date","time","strain","type"]);

G = groupsummary(data_tbl,["strain", "type"],"mean", ["magnitude_slow","magnitude_fast","sym_ratio_slow", "sym_ratio_fast", "DSI_vector_slow", "DSI_vector_fast", "DSI_pdnd_slow", "DSI_pdnd_fast"]);


%% 1 - Is there a difference in DSI between cells from different types / strains? 

data_type = 'sigma_y_exc';

%% Unbalanced n-way ANOVA:

[p, tbl, stats] = anovan(data_tbl.(data_type), {data_tbl.type, data_tbl.strain}, ...
                         'model', 'interaction', ...
                         'varnames', {'Type', 'Strain'});

c = multcompare(stats, 'Dimension', [1,2]);

%% Generate boxplot:

generate_boxplot_data_tbl(data_tbl, data_type)


%% Scatterplots
% Compare DSI with quality metrics

vals_ctl = find(data_tbl.strain == "control");
vals_ttl = find(data_tbl.strain ~= "control");

%%
data_type = "sym_ratio_fast";

figure
scatter(qual_res_table.var_v(vals_ctl), data_tbl.(data_type)(vals_ctl), 'ko');
hold on
scatter(qual_res_table.var_v(vals_ttl), data_tbl.(data_type)(vals_ttl), 'ro');
[rho, pval] = corr(qual_res_table.var_v, data_tbl.sym_ratio_fast, 'Type','Spearman')

xlabel('recording variance (mV)')
ylabel(data_type)
xlim([0 40])
title(strcat("pval = ", string(pval)))
ax = gca;
ax.TickDir = 'out';
f = gcf;
f.Position = [620   662   367   305];

%%
data_type = "sym_ratio_slow";

figure
scatter(qual_res_table.var_v(vals_ctl), data_tbl.(data_type)(vals_ctl), 'ko');
hold on
scatter(qual_res_table.var_v(vals_ttl), data_tbl.(data_type)(vals_ttl), 'ro');
[rho, pval] = corr(qual_res_table.var_v, data_tbl.(data_type), 'Type','Spearman')
xlabel('recording variance (mV)')
ylabel(data_type)
xlim([0 40])
title(strcat("pval = ", string(pval)))
ax = gca;
ax.TickDir = 'out';
f = gcf;
f.Position = [620   662   367   305];

%%
data_type = "DSI_vector_slow";

figure
scatter(qual_res_table.var_v(vals_ctl), data_tbl.(data_type)(vals_ctl), 'ko');
hold on
scatter(qual_res_table.var_v(vals_ttl), data_tbl.(data_type)(vals_ttl), 'ro');
[rho, pval] = corr(qual_res_table.var_v, data_tbl.(data_type), 'Type','Spearman')
xlabel('recording variance (mV)')
ylabel(data_type)
xlim([0 40])
title(strcat("pval = ", string(pval)))
ax = gca;
ax.TickDir = 'out';
f = gcf;
f.Position = [620   662   367   305];


%%
data_type = "DSI_vector_fast";

figure
scatter(qual_res_table.var_v(vals_ctl), data_tbl.(data_type)(vals_ctl), 'ko');
hold on
scatter(qual_res_table.var_v(vals_ttl), data_tbl.(data_type)(vals_ttl), 'ro');
[rho, pval] = corr(qual_res_table.var_v, data_tbl.(data_type), 'Type','Spearman')
xlabel('recording variance (mV)')
ylabel(data_type)
xlim([0 40])
title(strcat("pval = ", string(pval)))
ax = gca;
ax.TickDir = 'out';
f = gcf;
f.Position = [620   662   367   305];


%%
data_type = "DSI_vector_fast";
data_type_x = "R_squared";

figure
scatter(rf_res_table.(data_type_x)(vals_ctl), data_tbl.(data_type)(vals_ctl), 'ko');
hold on
scatter(rf_res_table.(data_type_x)(vals_ttl), data_tbl.(data_type)(vals_ttl), 'ro');
[rho, pval] = corr(rf_res_table.(data_type_x), data_tbl.(data_type), 'Type','Spearman')
xlabel(strrep(data_type_x, '_', '-'))
ylabel(strrep(data_type, '_', '-'))
% xlim([0 40])
title(strcat("pval = ", string(pval)))
ax = gca;
ax.TickDir = 'out';
f = gcf;
f.Position = [620   662   367   305];



%% QUALITY OF RECORDINGS

data_tbl = qual_res_table;
data_tbl = renamevars(data_tbl,["Date","Time","Strain","Type"],["date","time","strain","type"]);

G = groupsummary(data_tbl,["strain", "type"],"mean", ["median_v", "var_v", "max_v", "min_v"]);


%% Unbalanced n-way ANOVA:

data_type = "min_v";

[p, tbl, stats] = anovan(data_tbl.(data_type), {data_tbl.type, data_tbl.strain}, ...
                         'model', 'interaction', ...
                         'varnames', {'Type', 'Strain'});

c = multcompare(stats, 'Dimension', [1,2]);

%% Generate boxplot:

generate_boxplot_data_tbl(data_tbl, data_type)




%% RF ESTIMATE

data_tbl = rf_res_table;
data_tbl = renamevars(data_tbl,["Date","Strain","Type"],["date","strain","type"]);

G = groupsummary(data_tbl,["strain", "type"],"mean", ["med_var_X_reps", "med_var_W_reps", "med_diff_mean", "R_squared", "R_squaredi"]);


%% Unbalanced n-way ANOVA:

data_type = "R_squaredi";

[p, tbl, stats] = anovan(data_tbl.(data_type), {data_tbl.type, data_tbl.strain}, ...
                         'model', 'interaction', ...
                         'varnames', {'Type', 'Strain'});

c = multcompare(stats, 'Dimension', [1,2]);

%% Generate boxplot:

generate_boxplot_data_tbl(data_tbl, data_type)










