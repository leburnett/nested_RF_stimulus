% Comparisons of RNAi and MARCM data

%% Mean Vm

data1 = cell2mat(summaryTable.mean_vm);
data2 = cell2mat(summary_MARCM.mean_vm);
ttl = 'Mean Vm';


%% STD

data1 = cell2mat(summaryTable.std_vm);
data2 = cell2mat(summary_MARCM.std_vm);

ttl = 'STD Vm';

%% Max V

data1 = cell2mat(summaryTable.max_vm);
data2 = cell2mat(summary_MARCM.max_vm);
ttl = 'Max Vm';

%% Min V

data1 = cell2mat(summaryTable.min_vm);
data2 = cell2mat(summary_MARCM.min_vm);
ttl = 'Min Vm';


%% Control RNAi versus TTL RNAi %%  %%% % %%% %%% %% %% %%% %% 

%% Mean Vm 

rows_ctrl = strcmp(string(summaryTable.comments), "control RNAi");
data1 = cell2mat(summaryTable.mean_vm(rows_ctrl));

rows_ttl = strcmp(string(summaryTable.comments), 'turtle RNAi');
data2 = cell2mat(summaryTable.mean_vm(rows_ttl));

ttl = 'Mean Vm';


%% STD_vm

rows_ctrl = strcmp(string(summaryTable.comments), "control RNAi");
data1 = cell2mat(summaryTable.std_vm(rows_ctrl));

rows_ttl = strcmp(string(summaryTable.comments), 'turtle RNAi');
data2 = cell2mat(summaryTable.std_vm(rows_ttl));

ttl = 'STD Vm';


%% Max V

rows_ctrl = strcmp(string(summaryTable.comments), "control RNAi");
data1 = cell2mat(summaryTable.max_vm(rows_ctrl));

rows_ttl = strcmp(string(summaryTable.comments), 'turtle RNAi');
data2 = cell2mat(summaryTable.max_vm(rows_ttl));

ttl = 'Max Vm';

%% Min V

rows_ctrl = strcmp(string(summaryTable.comments), "control RNAi");
data1 = cell2mat(summaryTable.min_vm(rows_ctrl));

rows_ttl = strcmp(string(summaryTable.comments), 'turtle RNAi');
data2 = cell2mat(summaryTable.min_vm(rows_ttl));

ttl = 'Min Vm';


S = compareTwoGroupsNP(data1, data2,'labels', {'ctrl RNAi','ttl RNAi'}, 'paired', false, 'title', ttl, 'ylabel', 'Voltage (mV)', 'xlabel', '', 'colors', [0.5 0.5 0.5; 1 0.5 0.5]);

%% Control vs TTL in MARCM


%% Mean Vm 

rows_ctrl = contains(string(summary_MARCM.comments), "control");
data1 = cell2mat(summary_MARCM.mean_vm(rows_ctrl));

rows_ttl = contains(string(summary_MARCM.comments), "tut mut");
data2 = cell2mat(summary_MARCM.mean_vm(rows_ttl));

ttl = 'Mean Vm';


%% STD_vm

rows_ctrl = contains(string(summary_MARCM.comments), "control");
data1 = cell2mat(summary_MARCM.std_vm(rows_ctrl));

rows_ttl = contains(string(summary_MARCM.comments), "tut mut");
data2 = cell2mat(summary_MARCM.std_vm(rows_ttl));

ttl = 'STD Vm';


%% Max V

rows_ctrl = contains(string(summary_MARCM.comments), "control");
data1 = cell2mat(summary_MARCM.max_vm(rows_ctrl));

rows_ttl = contains(string(summary_MARCM.comments), "tut mut");
data2 = cell2mat(summary_MARCM.max_vm(rows_ttl));

ttl = 'Max Vm';

%% Min V

rows_ctrl = contains(string(summary_MARCM.comments), "control");
data1 = cell2mat(summary_MARCM.min_vm(rows_ctrl));

rows_ttl = contains(string(summary_MARCM.comments), "tut mut");
data2 = cell2mat(summary_MARCM.min_vm(rows_ttl));

ttl = 'Min Vm';


S = compareTwoGroupsNP(data1, data2,'labels', {'ctrl MARCM','ttl MARCM'}, 'paired', false, 'title', ttl, 'ylabel', 'Voltage (mV)', 'xlabel', '', 'colors', [0.5 0.5 0.5; 1 0.5 0.5]);






















