
T2 = combine_bar_results();

% Add Genotype column. 
T.Genotype{1}='control';


%% DSI vector

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = (T.DSI_vector(rows_ctrl));

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = (T.DSI_vector(rows_ttl));

ttl = 'DSI-vector';


%% DSI pdnd

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = (T.DSI_pdnd(rows_ctrl));

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = (T.DSI_pdnd(rows_ttl));

ttl = 'DSI-pdnd';


%% Median voltage

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.median_voltage(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.median_voltage(rows_ttl);

ttl = 'median-voltage';


%% Magnitude of vector sum

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.magnitude(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.magnitude(rows_ttl);

ttl = 'magnitude';


%%
rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.fwhm(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.fwhm(rows_ttl);

ttl = 'fwhm';


%% Symmetry ratio

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.sym_ratio(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.sym_ratio(rows_ttl);

ttl = 'sym-ratio';


%% CV

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.cv(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.cv(rows_ttl);

ttl = 'CV';


%% thetahat

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.thetahat(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.thetahat(rows_ttl);

ttl = 'thetahat';

%% kappa

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.kappa(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.kappa(rows_ttl);

ttl = 'kappa';


%% v_max

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.v_max(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.v_max(rows_ttl);

ttl = 'v-PD';
ylb = 'Voltage (mV)';


%% v_null

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.v_null(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.v_null(rows_ttl);

ttl = 'v-ND';
ylb = 'Voltage (mV)';



%% v_ortho1

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.v_ortho1(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.v_ortho1(rows_ttl);

ttl = 'v-ortho1';
ylb = 'Voltage (mV)';


%% v_ortho2

rows_ctrl = strcmp((T.Genotype), 'control');
data1 = T.v_ortho2(rows_ctrl);

rows_ttl = strcmp((T.Genotype), 'turtle');
data2 = T.v_ortho2(rows_ttl);

ttl = 'v-ortho2';
ylb = 'Voltage (mV)';



S = compareTwoGroupsNP(data1, data2,'labels', {'ctrl MARCM','ttl MARCM'}, 'paired', false, 'title', ttl, 'ylabel', ylb, 'xlabel', '', 'colors', [0.5 0.5 0.5; 1 0.5 0.5]);


























%% RNAi data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% DSI vector

rows_ctrl = strcmp((T.Strain), 'control');
data1 = cell2mat(T.DSI_vector(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = cell2mat(T.DSI_vector(rows_ttl));

ttl = 'DSI-vector';
ylb = 'DSI';


%% DSI pdnd

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.DSI_pdnd(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.DSI_pdnd(rows_ttl));

ttl = 'DSI-pdnd';
ylb = 'DSI';


%% Median voltage

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.median_voltage(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.median_voltage(rows_ttl));

ttl = 'median-voltage';
ylb = 'Voltage (mV)';


%% Magnitude of vector sum

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.magnitude(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.magnitude(rows_ttl));

ttl = 'magnitude';


%% Symmetry ratio

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.sym_ratio(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.sym_ratio(rows_ttl));

ttl = 'sym-ratio';
ylb = 'Ratio';


%% CV

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.cv(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.cv(rows_ttl));

ttl = 'CV';
ylb = 'Circular Variance';


%% thetahat

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.thetahat(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.thetahat(rows_ttl));

ttl = 'thetahat';

%% kappa

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.kappa(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.kappa(rows_ttl));

ttl = 'kappa';
ylb = 'Concentration Parameter';


%% v_max
rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.v_max(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.v_max(rows_ttl));

ttl = 'v-PD';
ylb = 'Voltage (mV)';

%% v_null

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.v_null(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.v_null(rows_ttl));

ttl = 'v-ND';
ylb = 'Voltage (mV)';


%% v_ortho1

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.v_ortho1(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.v_ortho1(rows_ttl));

ttl = 'v-ortho1';
ylb = 'Voltage (mV)';

%% v_ortho2

rows_ctrl = strcmp((T.Strain), 'control');
data1 = (T.v_ortho2(rows_ctrl));

rows_ttl = strcmp((T.Strain), 'ttl');
data2 = (T.v_ortho2(rows_ttl));

ttl = 'v-ortho2';
ylb = 'Voltage (mV)';


S = compareTwoGroupsNP(data1, data2,'labels', {'ctrl RNAi','ttl RNAi'}, 'paired', false, 'title', ttl, 'ylabel', ylb, 'xlabel', '', 'colors', [0.5 0.5 0.5; 1 0.5 0.5]);
