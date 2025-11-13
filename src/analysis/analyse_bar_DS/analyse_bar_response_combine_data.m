
% Assess bar responses across the different cohorts - split into ON / OFF. 
% Created Nov 5th 2025

%% 1 - Create the table:
T = combine_bar_results();
T = addOrthoMetrics(T);

T_on = T(T.Type == "on", :);
T_off = T(T.Type == "off", :);


%% 2 - Polar plots:

polar_mean_by_strain(T_on, "slow")
title("MARCM-ttl-control-ON-w+", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];


polar_mean_by_strain(T_off, "slow")
title("MARCM-ttl-control-OFF-w+", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];

%% 3 - Box and whisker plots:

columnName = "DSI_vector";
plotGroupedBox(T, columnName)

%% 4 - Stats

results = runGroupedStats(T, columnName);


%%  5 - Polar Scatter  - metrics versus PD angle.

plotPolarByGroup(T, "v_null_v_max_slow")
f = gcf;
f.Position = [620   621   457   346];

rlim([0 1])


%% Polar plots - shape based on PD angle. 

T_1 = T(T.angGroup == 1, :);
T_2 = T(T.angGroup == 2, :);
T_3 = T(T.angGroup == 3, :);
T_4 = T(T.angGroup == 4, :);

polar_mean_by_strain(T_4, "slow")
title("Group-4", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];

%% ON cells grouped by PD 

T_1 = T_on(T_on.angGroup == 1, :);
T_2 = T_on(T_on.angGroup == 2, :);
T_3 = T_on(T_on.angGroup == 3, :);
T_4 = T_on(T_on.angGroup == 4, :);

polar_mean_by_strain(T_1, "slow")
title("Group-1-ON", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];

%% OFF cells grouped by PD 

T_1 = T_off(T_off.angGroup == 1, :);
T_2 = T_off(T_off.angGroup == 2, :);
T_3 = T_off(T_off.angGroup == 3, :);
T_4 = T_off(T_off.angGroup == 4, :);

polar_mean_by_strain(T_4, "slow")
title("Group-4-OFF", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];






















































%%

% % For the RNAi data:
% % "251017_all_results_table_RNAi_ttl.mat"
% polar_mean_by_strain(T)
% title("RNAi-ttl-control-ON", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];
% 
% % For the MARCM data:
% % "251017_all_results_table_MARCM_ttl.mat"
% polar_mean_by_genotype(T)
% title("MARCM-ttl-control-OFF", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];


%% For the MARCM data from Autumn 2025 onwards. "control" / "ttl" now added in the strain name. 
% 
% polar_mean_by_strain(T, "slow")
% title("MARCM-ttl-control-OFF-w+", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];

% Convert from cell to double:

% % Check the variable types first
% class(T.DSI_vector)
% class(T_new.DSI_vector)
% 
% % If any of them is a cell array, convert it to a numeric array
% if iscell(T.DSI_vector)
%     T.DSI_vector = cell2mat(T.DSI_vector);
% end
% 
% if iscell(T_new.DSI_vector)
%     T_new.DSI_vector = cell2mat(T_new.DSI_vector);
% end
% 
% % Now combine the tables
% T_combined = [T; T_new];

% strainSel = "control";
% polar_mean_by_type_for_strain(T, strainSel)
% 
% strainSel = "ttl";
% polar_mean_by_type_for_strain(T, strainSel)


%% Different speeds

% 
% polar_mean_by_strain(T, "slow")
% title("MARCM-ttl-control-OFF-28dps", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];
% 
% polar_mean_by_strain(T, "fast")
% title("MARCM-ttl-control-ON-56dps", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];
% 
% polar_mean_by_strain(T, "vfast")
% title("MARCM-ttl-control-OFF-168dps", 'FontSize', 16)
% f = gcf;
% f.Position = [620   510   620   457];





























