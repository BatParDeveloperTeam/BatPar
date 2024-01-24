

%% mode selection
flags.chagre_mode.tf = evalin('base','Mode_charge_tf');
flags.cathode_mode.tf = evalin('base','Mode_cathode_tf'); % only used in easyECMfit for volatge constraints
flags.anode_mode.tf = evalin('base','Mode_anode_tf');

if strcmp(flags.chagre_mode.tf,'yes')
    if strcmp(flags.anode_mode.tf,'yes')
        working_mode = 1;   % mode 1: anode - charge
    else
        working_mode = 2;   % mode 2: full cell/cathode - charge
    end 
else 
    if strcmp(flags.anode_mode.tf,'yes')
        working_mode = 3;   % mode 3: anode - discharge
    else 
        working_mode = 4;   % mode 4: full cell/cathode - discharge
    end 
end 

%% More switches

flags.R0temp.tf   = 'no';     % do you want R0 to depend on temperature? Please select no, as parameterisation tests are normally perfromed under invariable temperature (ideally)

NRC = str2num(evalin('base','RCpair_number'));                      % how many R-C pairs do you want (at least 1)

flags.taufixed_Stage2.tf = evalin('base','fixedtau_realy_tf');     % Will there be a Stage 2 for solving fixed tau

flags.taufix.tf  = 'no'; % (this indicator only works in Stage 1) Initially, tau is not fixed in the frist stage.

flags.pulsehead4R0.tf = evalin('base','Pulse_head_tf');

flags.allplots.tf = 'no';      % do you want to see individual fits (with optimised initial voltages over RC pairs)?

flags.Ritaui_relimit.tf  = evalin('base','Ri_taui_relimit_tf'); % do you want to relimit the upper limits of tai_i and R_i
% flags.deeper_datasplitter.tf  = evalin('base','deeperOP_tf'); % do you want to squeeze the range of data split to each SOC

flags.OCVpseudo.tf = evalin('base','OCVpseudo_tf'); % Parameterisation of pseudo OCV
flags.OCVtrue.tf = evalin('base','OCVtrue_tf'); % Parameterisation of true OCV
flags.R0only.tf = evalin('base','R0only_tf'); % Parameterisation of true OCV

% current and/or temperature dependence
flags.AmpDependence.tf = evalin('base','CurrentDenp_tf');              % current dependence
Currents_depended = sort(str2num(evalin('base','Currents_4_Denp')))';       % currents (absolute value)
Currents_sensitivity = str2num(evalin('base','CurrentSensi_4_Denp'));  % current sensitivity

flags.TempDependence.tf = evalin('base','TempDenp_tf');                % Temperature dependence
Temps_depended = sort(str2num(evalin('base','Temps_4_Denp')))';               % Temperatures
Temps_sensitivity = str2num(evalin('base','TempSensi_4_Denp'));        % Temperature sensitivity

%% Initialisation - read and sort data 

cd(evalin('base','exp_path'))               % goes into the data folder

M = xlsread(evalin('base','datasheet_name'));
M = M (~isnan(M(:,1)),:);
N = M(:,[str2num(evalin('base','time_column')),str2num(evalin('base','current_column')),str2num(evalin('base','voltage_column')),str2num(evalin('base','record_column')),str2num(evalin('base','temp_column'))]); % keep 5 columns - time, current, voltage, Rec#, temperature
[~,ind,~] = unique(N(:,1)); % ~:neglected  ind:row mubers of unique times in N
N = N(ind,:);               % keep data at unique times

flags.datasheet_cleanup.tf = evalin('base','datasheet_cleanup_tf');
if strcmp(flags.datasheet_cleanup.tf,'yes')
exp_start_line = str2num(evalin('base','row_start')); % the line number in N, where the experiment formally strats (no discharge or chagre preparations)
exp_end_line = str2num(evalin('base','row_end')); % the line number in N, where the experiment formally ends (no discharge or chagre preparations)   

    if isempty(exp_start_line)
        exp_start_line = 1;
    else
        exp_start_line = find (N(:,4) == exp_start_line);
    end 
    if isempty(exp_end_line)
        [rowsN,columnsN]=size (N);
        exp_end_line = rowsN;
    else 
        exp_end_line = find (N(:,4) == exp_end_line);
    end
else
    exp_start_line = 1;
    [rowsN,columnsN]=size (N);
    exp_end_line = rowsN;
end 

flags.dischargecurrentsubzero.tf = evalin('base','Current_belwozero_tf');
if strcmp(flags.dischargecurrentsubzero.tf,'yes')
    Current_sign = -1;
else
    Current_sign = 1;
end

t_ = N(exp_start_line:exp_end_line,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_ = Current_sign * N(exp_start_line:exp_end_line,2)/1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if working_mode == 1 || working_mode == 3 % anode paramterisation
    anode_voltage_sign = mean(N(exp_start_line:exp_end_line,3));
    if anode_voltage_sign < 0 % input andoe voltage is mostly below zero
    V_ = N(exp_start_line:exp_end_line,3);
    else  % input andoe voltage is mostly above zero
    V_ = - N(exp_start_line:exp_end_line,3); 
    end 
else
    V_ = N(exp_start_line:exp_end_line,3);
end
T_ = N(exp_start_line:exp_end_line,5);

cd ..
%% SOC calculation and processing

Ccap = str2num(evalin('base','Nom_capacity'))*3600;                % 4.4 Amp hours
S_ = - cumtrapz(t_,I_/Ccap);  % delta SOC from time zero; calculated by coulomb counting

bounds.upperbound = nan;  % the struct to save upper and lower bounds of SOC; to be used in datasplitter.m
bounds.lowerbound = nan;

flags.MaxSOCis1.tf = evalin('base','MaxSOCis1'); % 100% SOC is the baseline
if strcmp(flags.MaxSOCis1.tf,'yes')
    S_real = S_ + (1-max(S_)); % Real SOC at each time; (1-max(S_))=initial SOC;
    S_= S_real;
    % S_(find(S_<0)) = 0; % Real SOC filters out minus values
    S_lowerbound = min(S_); % minimum SOC to be parameterised
    bounds.lowerbound = S_lowerbound;
else  % 0% SOC is the baseline
    S_real = S_ + (0-min(S_)); % Real SOC at each time; (0-min(S_))=initial SOC;
    S_= S_real;
    % S_(find(S_>1)) = 1; % Real SOC filters out values that outnumber 1
    S_upperbound = max(S_); % minimum SOC appeared in experiment
    bounds.upperbound = S_upperbound;
end

%% optional - import SOC - OCV table 
flags.SOCOCVtable.tf   = evalin('base','OCVtable_tf');     % do you have an experimental SOC-OCV table to import ?
% The table should be prepared as: two columns, the first one being SOC,the second one being OCV

SOCTab= [];
OCVTab = [];
if (strcmp(flags.SOCOCVtable.tf,'yes'))  % string compare
cd SOC_OCV_table    
P = xlsread(evalin('base','OCVsheet_name'));
Q = P; %
[~,ind_,~] = unique(Q(:,1)); % ~:neglected  ind:row mubers of unique times in Q
Q = Q(ind_,:);               % keep data at unique SOC

SOCTab = Q(:,1);
OCVTab = Q(:,2);
cd ..
end

%% optional - import SOC - R0 table 
flags.SOCR0table.tf   = evalin('base','R0table_tf');     % do you have an experimental SOC-R0 table to import ?
% The table should be prepared as: two columns, the first one being SOC,the second one being R0

SOCSOCTab= [];
R0R0Tab = [];
if (strcmp(flags.SOCR0table.tf,'yes'))  % string compare
cd SOC_R0_table    
PP = xlsread(evalin('base','R0sheet_name'));
QQ = PP; %
[~,ind_,~] = unique(QQ(:,1)); % ~:neglected  ind:row mubers of unique times in Q
QQ = QQ(ind_,:);               % keep data at unique SOC

SOCSOCTab = QQ(:,1);
R0R0Tab = QQ(:,2);
cd ..
end

%% discharge/chagre data sortation
[t, I, V, S, T, t_charge, I_charge, V_charge, S_charge, T_charge] = dataextractor (t_, I_, V_, S_, T_); % Extract continuous discharge (and charge) segments from experiment data

%% OCV or R0 ONLY
if strcmp(flags.R0only.tf,'yes')
    disp('Parameterisation of R0 ONLY')    
end

if strcmp(flags.OCVtrue.tf,'yes')
    disp('Parameterisation of true OCV ONLY')
    disp('Please be noted - ture OCV parameterisation takes into account both DISCHAGRE and CHARGE data, if you have both in your datasheet') 
end

if strcmp(flags.OCVpseudo.tf,'yes')
    disp('Parameterisation of pseudo OCV ONLY') 
    disp('Please be noted - pseudo OCV parameterisation takes into account either DISCHAGRE or CHARGE data, based on your selection in the GUI (or in {list.xlsx})') 
end

%% plot discharge & charge segments
if strcmp(flags.OCVtrue.tf,'no')
    if working_mode == 3 || working_mode == 4
    disp([num2str(numel(t)),' continuous segment(s) of DISCHARGE extracted from your datasheet.'])
    else
    disp([num2str(numel(t_charge)),' continuous segment(s) of CHARGE extracted from your datasheet.'])
    end 
end

cd Results
cd(fileNAME)

if working_mode == 3 || working_mode == 4
f1=figure;
hold on
for kk = 1:numel(t) % discharge segments
    plot(S{kk,1},V{kk,1})
    xlabel('SOC');
    ylabel('Voltage (V)');
    title (['Data extraction - Voltage VS SOC of ',num2str(numel(t)),' discharge segments'])
end
savefig('Data extraction of discharge segments.fig')

else
f1=figure;
hold on
for kk = 1:numel(t_charge) % charge segments
    plot(S_charge{kk,1},V_charge{kk,1}) 
    xlabel('SOC');
    ylabel('Voltage (V)');
    title (['Data extraction - Voltage VS SOC of ',num2str(numel(t_charge)),' charge segments'])
end
savefig('Data extraction of charge segments.fig') 

end
cd ..
cd ..
