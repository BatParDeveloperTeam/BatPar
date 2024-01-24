clear all
clearvars

warning('off', 'all');

tic % timmer starts

fileNAME = evalin('base','datasheet_name');

% naming of the result folder(s)
cathodeORnot = evalin('base','Mode_cathode_tf');
anodeORnot = evalin('base','Mode_anode_tf');
if strcmp(cathodeORnot,'yes')                
    fileNAME = [fileNAME,'_Cathode'];
else
    if strcmp(anodeORnot,'yes')  
        fileNAME = [fileNAME,'_Anode'];
    else
        fileNAME = [fileNAME,'_FullCell'];
    end
end

chargeORnot = evalin('base','Mode_charge_tf');
if strcmp(chargeORnot,'yes')
    fileNAME = [fileNAME,'_CHARGE'];
else
    fileNAME = [fileNAME,'_DISCHARGE'];
end

currentORnot= evalin('base','CurrentDenp_tf');
tempORnot= evalin('base','TempDenp_tf');
if strcmp(currentORnot,'yes') &&  strcmp(tempORnot,'no')
    fileNAME = [fileNAME,'_CrateSplit'];
end
if strcmp(currentORnot,'no') &&  strcmp(tempORnot,'yes')
    fileNAME = [fileNAME,'_TempSplit'];
end
if strcmp(currentORnot,'yes') &&  strcmp(tempORnot,'yes')
    fileNAME = [fileNAME,'_CrateTempSplit'];
end

tureOCVORnot= evalin('base','OCVtrue_tf');
pseudoOCVORnot= evalin('base','OCVpseudo_tf');
R0ORnot= evalin('base','R0only_tf');

if strcmp(tureOCVORnot,'yes')
    fileNAME = [fileNAME,'_tureOCVonly'];
end
if strcmp(pseudoOCVORnot,'yes') 
    fileNAME = [fileNAME,'_pseudoOCVonly'];
end
if strcmp(R0ORnot,'yes') 
    fileNAME = [fileNAME,'_R0only'];
end


cd Results
mkdir(fileNAME)
cd(fileNAME)
diary 'Command window log.txt' % log file starts
cd ..
cd ..

disp(['Directing to ', evalin('base','datasheet_name')])
disp('Initialising...')
%% Load in experimental data
 
 loadexperiments_

%% Parameterisation of OCV ONLY

if strcmp(flags.OCVtrue.tf,'yes') || strcmp(flags.OCVpseudo.tf,'yes')
    if strcmp(flags.AmpDependence.tf,'yes') || strcmp(flags.TempDependence.tf,'yes')
        disp('Please disable Current and/or Temperature dependence(s), since OCV-ONLY parameterisation along with those dependences has not been developed')
        return
    else 
        OCVparameterisationONLY
    end 
else
%% Division of SOC windows

S_interval = str2num(evalin('base','SOC_window'));    % the SOC interval you want in look-up tables;  >= 0.01

if strcmp(flags.MaxSOCis1.tf,'yes')  % 100% SOC is the baseline
S_lowerlimit = ceil (S_lowerbound/S_interval)*S_interval; % the lower limit of DIVIDABLE SOC in look-up tables
S_j = S_lowerlimit:S_interval:1; % SOC points: from S_lowerlimit to 1, with a fixed interval - S_interval
S_j = (vertcat(S_lowerbound,S_j'))'; % SOC points: from S_lowerbound to S_lowerlimit to 1
[~,ind,~] = unique(S_j(1,:)); % ~:neglected  ind:row mubers of unique times in N
S_j = S_j(:,ind);               % keep data at unique times
else   % 0% SOC is the baseline
S_upperlimit = 1- ceil ((1- S_upperbound)/S_interval)*S_interval; % the upper limit of DIVIDABLE SOC in look-up tables
S_j = 0:S_interval:S_upperlimit; % SOC points: from 0 to S_upperlimit, with a fixed interval - S_interval
S_j = (vertcat(S_j', S_upperbound))'; % SOC points: from 0 to S_upperlimit to S_upperbound
[~,ind,~] = unique(S_j(1,:)); % ~:neglected  ind:row mubers of unique times in N
S_j = S_j(:,ind);
end 

% Specify the lower and upper limits for SOC
SOC_param_lowlimit = str2num(evalin('base','SOC_param_lowlimit'));
if isempty(SOC_param_lowlimit) == 0
    S_j = S_j(find(S_j > SOC_param_lowlimit));
    S_j = (vertcat(SOC_param_lowlimit,S_j'))';
    bounds.lowerbound = SOC_param_lowlimit;
end

SOC_param_uplimit = str2num(evalin('base','SOC_param_uplimit'));
if isempty(SOC_param_uplimit) == 0
    S_j = S_j(find(S_j < SOC_param_uplimit));
    S_j = (vertcat(S_j',SOC_param_uplimit))';
    bounds.upperbound = SOC_param_uplimit;
end
%%%%

S_j = roundn(S_j,-6);  % update on 18 Apr 2022 --- to dismiss the very tiny difference (e.g. 1e-17) caused by 'str2num'
S_j = unique(S_j);     % update on 18 Apr 2022 --- to dismiss the very tiny difference (e.g. 1e-17) caused by 'str2num'

%% Preporcess - data sorted into each SOC window and clean up
if working_mode == 3 || working_mode == 4
[data, S_j] = datasplitter(t,I,V,T,S,S_j,bounds,flags); % sort out experimental data as SOC dependent
else

[data, S_j] = datasplitter(t_charge,I_charge,V_charge,T_charge,S_charge,S_j,bounds,flags); % sort out experimental data as SOC dependent
end

disp(['SOC range considered in parameterisation is between ',num2str(min(S_j)),' and ', num2str(max(S_j)),'.'])

if strcmp(flags.pulsehead4R0.tf,'yes')
    disp('Pulse Head (rather than Pulse End) is used to parameterise R0.')
else
    disp('Pulse End (rather than Pulse Head) is used to parameterise R0.')
end

if strcmp(flags.SOCOCVtable.tf,'yes')
    disp('SOC-OCV table is prescribed rather being parameterised.')
end

if strcmp(flags.SOCR0table.tf,'yes')
    disp('SOC-R0 table is prescribed rather being parameterised.')
end

% remove very short sets of data 
data = cleandown(data,NRC); 

%%% Belwo not useful commands, but leave them in case of unexpected incompatibilities
flags = yesnoconvert(flags); % yes & no converted into 1 & 0

flags.R0temp.val = [0,100];  % read off temperatures for R0

%% Switch - If current and/or temperature dependence(s) needed, codes will be redirected
if flags.AmpDependence.tf == 1 && flags.TempDependence.tf == 0
    AmpDenpOnly
else
    if flags.AmpDependence.tf == 0 && flags.TempDependence.tf == 1
        TempDenpOnly
    else    
        if flags.AmpDependence.tf == 1 && flags.TempDependence.tf == 1
            AmpAndTempDenp
        else
%% From here down, the codes only work for the case where neither current nor temperature dependence is considered
if flags.R0only.tf == 0
    if flags.taufixed_Stage2.tf == 0
    disp(['Initialisation done - Paramterisation (with VARIABLE tau, without I or T dependence) in progress, including ',num2str(numel(S_j)),' STEPS.']) 
    disp('***************************************************************************************************')
    else
    disp(['Initialisation done - Paramterisation (with constant tau, without I or T dependence) in progress, including 2 STAGES --> Each STAGE further includes ',num2str(numel(S_j)),' STEPS.']) 
    disp('***************************************************************************************************')
    disp('--> Proceeding to Stage 1 - solving parameters with VARIABLE tau.') 
    end
else
    disp(['Initialisation done - Paramterisation (without I or T dependence) in progress, including ',num2str(numel(S_j)),' STEPS.']) 
    disp('***************************************************************************************************')
end
%% Iterate over LUT - Stage 1 - solving wtih variable tau

% Initialisation of unkonws (variables to be paramterised)
OCV = zeros(1,numel(S_j))*nan;   
R0  = zeros(1+flags.R0temp.tf,numel(S_j))*nan;
R   = zeros(NRC,numel(S_j))*nan;
tau = zeros(NRC,numel(S_j))*nan;

%---------- Initialisation for the calculation of RMSE
L_temp = []; % wipe out values assigned in previous operation
for tao = 1:numel(S_j)
L_temp(tao) = sum(cellfun(@length, data.t{tao})); % sum length of all martrixes of one cell of data.t
end
L_dim = max(L_temp); % figure out the max(sum length of all martrixes of each cell of data.t)

DIFF = zeros(L_dim,numel(S_j))*nan; % the matrix records volatge difference between experiment and parameterised model

RMSE = zeros(1,numel(S_j)); % the vector records the root-mean-square error between experiment and model at different S-j
%----------

%%% The show starts
for ii = 1:numel(S_j)
    
    disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j)),', parameterising at SOC ', num2str(S_j(ii)),])  % display the progress of parameterisation
                       
    t_ = data.t{ii}; 
    I_ = data.I{ii}; 
    V_ = data.V{ii};  % pass data
    T_ = data.T{ii}; 
    S_ = data.S{ii}; 
    
    [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
    % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
    
    if Skip == 1
        DIFF(1,ii) = nan; % record Diff under S_j(ii)
        RMSE(1,ii)= nan; % record rmse under S_j(ii)
        
        disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
        disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')

    else
        
        DIFF(1:length(Diff),ii) = Diff; % record Diff under S_j(ii)
        RMSE(1,ii)= rms(Diff); % record rmse under S_j(ii)
        if flags.R0only.tf == 0
            disp(['     Step ',num2str(ii),' done. Root-mean-square error of voltage is ',num2str(RMSE(1,ii)),])
        end
    end
end

%% Recommend fixed tau constants

tau_1_reco = prctile (tau(1,:),50);
if NRC >= 2
tau_2_reco = prctile (tau(2,:),50);  
else
tau_2_reco = [];   
end
if NRC >= 3
tau_3_reco = prctile (tau(3,:),50);
else 
tau_3_reco = [];
end

tau_recommended =[tau_1_reco, tau_2_reco,tau_3_reco]; % recommend fixed tau for rerunning the codes

%% Suggestion for tau & R relimitation

    RandTauRelimit

%% Stage 2 - solving fixed tau

if flags.taufixed_Stage2.tf == 1 && flags.R0only.tf == 0
    disp('--> Proceeding to Stage 2 - solving parameters with CONSTANT tau') 
     
    flags.taufix.tf  = 1;  % Enabling fixed tau functions
    flags.taufix.val  = tau_recommended; % Fixed tau from Stage 1
    
    for ii = 1:numel(S_j)
        
        disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j)),', parameterising at SOC ', num2str(S_j(ii)),])  % display the progress of parameterisation
                           
        t_ = data.t{ii}; 
        I_ = data.I{ii}; 
        V_ = data.V{ii};  % pass data
        T_ = data.T{ii}; 
        S_ = data.S{ii}; 
        flags.R0solved.val  = R0(ii); % Solved R0 from Stage 1

        [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
        % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
        if Skip == 1
            DIFF(1,ii) = nan; % record Diff under S_j(ii)
            RMSE(1,ii)= nan; % record rmse under S_j(ii)
            
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
    
        else
            
            DIFF(1:length(Diff),ii) = Diff; % record Diff under S_j(ii)
            RMSE(1,ii)= rms(Diff); % record rmse under S_j(ii)
            
            disp(['     Step ',num2str(ii),' done. Root-mean-square error of voltage is ',num2str(RMSE(1,ii)),])
        end
    end

end

%%
Ci = tau ./ R; % calculate Ci



%% OCV and R0 check 

%%% monotonicity detection for OCV
diffOCV = diff(OCV);
indexOCV = find(diffOCV < 0);
if isempty(indexOCV) == 0
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
    disp(['     Parameterised OCV is not monotonic at SOC = ',num2str(S_j(indexOCV+1)), '    --> You MUST consider smoothing OCV.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%%% outlier detection for R0
Outlier_ind = isoutlier(R0,'movmedian',3);  % Hampel filter --> necessary but not sufficient
if any (Outlier_ind)
    Outlier_SOC = S_j(Outlier_ind);
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
    disp(['     Possible outlier(s) detected for R0 at SOC = ',num2str(S_j(Outlier_ind)), '    --> You MAY consider smoothing R0.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%% Final plots

cd Results
cd(fileNAME)

% Plots of OCV, R0, Ri, Ci, tau_i with SOC
f2=figure; plot (S_j, OCV, 'LineWidth',3); hold on; scatter(S_j, OCV,150,'k.'); xlabel('SOC'); ylabel('OCV (V)'); title ('Parameterised OCV with SOC'); savefig('Parameterised OCV.fig'); 
f3=figure; plot (S_j, R0, 'LineWidth',3);  hold on; scatter(S_j, R0,150,'k.'); xlabel('SOC'); ylabel('R0 (Ohm)'); title ('Parameterised R0 with SOC'); savefig('Parameterised R0.fig'); 
f4=figure; plot (S_j, R, 'LineWidth',3);  xlabel('SOC'); ylabel('Ri (Ohm)'); title ('Parameterised Ri with SOC'); legend ('R1 ','R2 ','R3 '); hold on; scatter(S_j, R,150,'k.','HandleVisibility','off'); savefig('Parameterised Ri.fig');
f5=figure; plot (S_j, Ci, 'LineWidth',3);  xlabel('SOC'); ylabel('Ci (F)'); title ('Parameterised Ci with SOC'); legend ('C1 ','C2 ','C3 '); hold on; scatter(S_j, Ci,150,'k.','HandleVisibility','off');savefig('Parameterised Ci.fig');
f6=figure; plot (S_j, tau, 'LineWidth',3); xlabel('SOC'); ylabel('tau_i'); title ('Parameterised tau_i with SOC'); legend ('tau1 ','tau2 ','tau3 '); hold on; scatter(S_j, tau,150,'k.','HandleVisibility','off'); savefig('Parameterised tau_i.fig');

cd ..
cd ..

    % Overall horizon - Experiment VS Paramterised Model
    DIFF_global = [];
    
    f7=figure; hold on; 
    
    if working_mode == 1 || working_mode == 3 % anode
        if anode_voltage_sign < 0
        plot(N(:,1), N(:,3), 'LineWidth',1); 
        else
        plot(N(:,1), -N(:,3), 'LineWidth',1); 
        end 
    else
        plot(N(:,1), N(:,3), 'LineWidth',1);
    end
    
    xlabel('Time (s)'); ylabel('Voltage (V)'); title ('Experiment VS Parameterised Model - Global horizon'); 

    if working_mode == 3 || working_mode == 4 % discharge  
      for kk = 1:numel(t)
        [Vhat, t_useful] = fullsolve(t{kk},I{kk},T{kk},S{kk},S_j,OCV,R0,R,tau,flags,NRC);
        if isnan(Vhat)
            DIFF_global{kk,1} = nan;
        else
        aaa= find (N(:,1) == t_useful(1)); bbb= find (N(:,1) == t_useful(end));
        if working_mode == 3 % anode
            if anode_voltage_sign < 0
            DIFF_global{kk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
            else
            DIFF_global{kk,1} = Vhat + N(aaa:bbb,3);  % Global difference of segment kk
            end
        else
        DIFF_global{kk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
        end
        plot(t_useful,Vhat,'LineWidth',1); % Parameterised model data
        end
      end
      for i = 2:numel(t)+1
        leg_str{i} = ['Model result ',num2str(i-1)];  
      end
    else  % charge
      for kk = 1:numel(t_charge)
        [Vhat, t_useful] = fullsolve(t_charge{kk},I_charge{kk},T_charge{kk},S_charge{kk},S_j,OCV,R0,R,tau,flags,NRC);
        if isnan(Vhat)
            DIFF_global{kk,1} = nan;
        else
        aaa= find (N(:,1) == t_useful(1)); bbb= find (N(:,1) == t_useful(end));
        if working_mode == 1 % anode
            if anode_voltage_sign < 0
            DIFF_global{kk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
            else
            DIFF_global{kk,1} = Vhat + N(aaa:bbb,3);  % Global difference of segment kk  
            end 
        else
        DIFF_global{kk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
        end
        plot(t_useful,Vhat,'LineWidth',1); % Parameterised model data
        end
      end
      for i = 2:numel(t_charge)+1
        leg_str{i} = ['Model result ',num2str(i-1)];  
      end
    end
    leg_str{1}='Experiment'; 
    legend(leg_str); % legends for model results


cd Results
cd(fileNAME)


savefig('Global comparison - Experiment VS Paramterised Model.fig');

%% local and global RMSE 

%-----Finalisation of local RMSE calculation
DIFF_vector = reshape(DIFF,[],1); % reshape DIFF matrix to a coclumn vector
DIFF_vector = DIFF_vector(~isnan(DIFF_vector)); % filter out NaN  

%------global RMSE
RMSE_Gtemp = cell2mat(DIFF_global);
RMSE_Gtemp = RMSE_Gtemp(~isnan(RMSE_Gtemp));
RMSE_GLOBAL= rms(RMSE_Gtemp); 

%% Save data and export LUT - Time to call it a day

% generate a look-up table including every parameter, in the format of .xlsx
LookUpTable_AllInOne = zeros(length(S_j),12)*nan;
LookUpTable_AllInOne (:,1) = S_j';
LookUpTable_AllInOne (:,2) = OCV';
LookUpTable_AllInOne (:,3) = R0';
[good,~] = size(R); [luck,~] = size(tau);
LookUpTable_AllInOne (:,4:3+good) = R';
LookUpTable_AllInOne (:,7:6+luck) = Ci';
LookUpTable_AllInOne (:,10:9+luck) = tau';

LookUpTable_AllInOne = mat2cell(LookUpTable_AllInOne,ones(1,length(S_j)),[1 1 1 1 1 1 1 1 1 1 1 1]);
title = {'SOC','OCV','R0','R1','R2','R3','C1','C2','C3','tau1','tau2','tau3'};
LookUpTable_AllInOne = [title;LookUpTable_AllInOne];
writecell(LookUpTable_AllInOne,'LookUpTable_AllInOne.xlsx');

if flags.R0only.tf == 0
    disp(['All done. The local and global root-mean-square errors (with respect to voltage) are ',num2str(rms(DIFF_vector)),' and ',num2str(RMSE_GLOBAL),'.'])
end
disp(['Now you can find parameterisation results in the folder \Results\',fileNAME])
disp('***************************************************************************************************')
% figures poped up before
clearvars -except working_mode Ccap data DIFF DIFF_vector DIFF_global RMSE_GLOBAL flags I I_ I_charge N NRC OCV R Ci R0 RMSE S_real S S_ S_charge  S_interval S_j bounds t T t_ T_ t_charge T_charge tau tau_recommended V V_ V_charge datasheet_name fixedtau_realy_tf Nom_capacity OCVsheet_name OCVtable_tf RCpair_number SOC_window voltage_up_limit voltage_low_limit fileNAME Currents_depended Currents_sensitivity dataIII S_j_I 
evalin('base', 'save(''INPUT.mat'')');
save OUTPUT

toc
diary off
cd ..
cd ..
        end 
    end 
end 
end
