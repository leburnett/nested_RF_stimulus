
% Assess bar responses across the different cohorts - split into ON / OFF. 
% Created Nov 5th 2025

T_all = T;

% For the RNAi data:
% "251017_all_results_table_RNAi_ttl.mat"
polar_mean_by_strain(T)
title("RNAi-ttl-control-ON", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];

% For the MARCM data:
% "251017_all_results_table_MARCM_ttl.mat"
polar_mean_by_genotype(T)
title("MARCM-ttl-control-OFF", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];


%% For the MARCM data from Autumn 2025 onwards. "control" / "ttl" now added in the strain name. 

polar_mean_by_strain(T_new)
title("MARCM-ttl-control-OFF", 'FontSize', 16)
f = gcf;
f.Position = [620   510   620   457];


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

strainSel = "control";
polar_mean_by_type_for_strain(T, strainSel)

strainSel = "ttl";
polar_mean_by_type_for_strain(T, strainSel)





























