

[dataTTT, SOC_steps] = data_temperatureDependent(data, Temps_depended, Temps_sensitivity, S_j, NRC);

if flags.R0only.tf == 0
    if flags.taufixed_Stage2.tf == 0
    disp(['Initialisation done - Paramterisation (with VARIABLE tau and T dependence) in progress, including ',num2str(length(Temps_depended)),' TEMPERATURES --> Each TEMPERATURE furhter includes ~',num2str(max(SOC_steps)),' STEPS.']) 
    else
    disp(['Initialisation done - Paramterisation (with CONSTANT tau and T dependence) in progress, including ',num2str(length(Temps_depended)),' TEMPERATURES --> Each TEMPERATURE further includes 2 STAGES --> Each STAGE further includes ~',num2str(max(SOC_steps)),' STEPS.']) 
    end
    disp('***************************************************************************************************')
else
    disp(['Initialisation done - Paramterisation (with T dependence) in progress, including ',num2str(length(Temps_depended)),' TEMPERATURES --> Each TEMPERATURE furhter includes ~',num2str(max(SOC_steps)),' STEPS.']) 
    disp('***************************************************************************************************')
end 

for iiii = 1:length(Temps_depended)

disp(['-> Parameterisation taking place at TEMPERATURE ',num2str(iiii),' of ',num2str(length(Temps_depended)),' -- ', num2str(Temps_depended(iiii)),'degC.']) 
  
if working_mode == 3 || working_mode == 4
    f111=figure;
    hold on
    for kk = 1:numel(T) % discharge segment
        indexTTT = ((Temps_depended(iiii) - Temps_sensitivity) <= T{kk}) + (T{kk} <= (Temps_depended(iiii) + Temps_sensitivity)) == 2; 
        plot(S{kk}(indexTTT),V{kk}(indexTTT))
        xlabel('SOC');
        ylabel('Voltage (V)');
        title (['Data extraction (discharge) - Voltage VS SOC @ ',num2str(Temps_depended(iiii)),'degC']);
    end
    
cd Results
cd(fileNAME)
newfile = [fileNAME,'_at_',num2str(Temps_depended(iiii)),'degC'];
mkdir(newfile)
cd(newfile)
    
savefig(['Data extraction of discharge segments @ ',num2str(Temps_depended(iiii)),'degC.fig'])

else
    f222=figure;
    hold on
    for kk = 1:numel(T_charge) % discharge segment
        indexTTT = ((Temps_depended(iiii) - Temps_sensitivity) <= T_charge{kk}) + (T_charge{kk} <= (Temps_depended(iiii) + Temps_sensitivity)) == 2;
        plot(S_charge{kk}(indexTTT),V_charge{kk}(indexTTT))
        xlabel('SOC');
        ylabel('Voltage (V)');
        title (['Data extraction (charge) - Voltage VS SOC @ ', num2str(Temps_depended(iiii)),'degC']);
    end

cd Results
cd(fileNAME)
newfile = [fileNAME,'_at_',num2str(Temps_depended(iiii)),'degC'];
mkdir(newfile)
cd(newfile)  

savefig(['Data extraction of charge segments @ ',num2str(Temps_depended(iiii)),'degC.fig'])
end


cd ..
cd ..
cd ..


if flags.taufixed_Stage2.tf == 1 && flags.R0only.tf == 0
disp('--> Proceeding to Stage 1 - solving parameters with VARIABLE tau.') 
end

S_j_I = dataTTT.S_j_I{iiii, 1}{1,1};  % the real SOC at the current Temperature

% Initialisation of unkonws (variables to be paramterised)
OCV = zeros(1,numel(S_j_I));   
R0  = zeros(1+flags.R0temp.tf,numel(S_j_I));
R   = zeros(NRC,numel(S_j_I));
tau = zeros(NRC,numel(S_j_I));

%---------- Initialisation for the calculation of RMSE
L_temp = [];   % wipe out values assigned in previous operation
for tao = 1:numel(S_j_I)
L_temp(tao) = sum(cellfun(@length, dataTTT.t{iiii,1}{1,1}{tao})); % sum length of all martrixes of one cell of data.t
end
L_dim = max(L_temp); % figure out the max(sum length of all martrixes of each cell of data.t)

DIFF = zeros(L_dim,numel(S_j_I))*nan; % the matrix records volatge difference between experiment and parameterised model

RMSE = zeros(1,numel(S_j_I)); % the vector records the root-mean-square error between experiment and model at different S-j
%----------

%%% The show starts
flags.taufix.tf  = 0;  % provisionally disabling fixed tau functions
for ii = 1:numel(S_j_I)
    
    disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j_I)),', parameterising at SOC ', num2str(S_j_I(ii)),])  % display the progress of parameterisation
                       
    t_ = dataTTT.t{iiii,1}{1,1}{ii};  
    I_ = dataTTT.I{iiii,1}{1,1}{ii}; 
    V_ = dataTTT.V{iiii,1}{1,1}{ii};  % pass data
    T_ = dataTTT.T{iiii,1}{1,1}{ii};
    S_ = dataTTT.S{iiii,1}{1,1}{ii};
    
    [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j_I(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
    % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
    
    if Skip == 1
        DIFF(1,ii) = nan; % record Diff under S_j(ii)
        RMSE(1,ii)= nan; % record rmse under S_j(ii)
        
        disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j_I(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
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
    
    for ii = 1:numel(S_j_I)
        
        disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j_I)),', parameterising at SOC ', num2str(S_j_I(ii)),])  % display the progress of parameterisation
                           
         t_ = dataTTT.t{iiii,1}{1,1}{ii};  
         I_ = dataTTT.I{iiii,1}{1,1}{ii}; 
         V_ = dataTTT.V{iiii,1}{1,1}{ii};  % pass data
         T_ = dataTTT.T{iiii,1}{1,1}{ii};
         S_ = dataTTT.S{iiii,1}{1,1}{ii};
         flags.R0solved.val  = R0(ii); % Solved R0 from Stage 1

        [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j_I(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
        % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
        if Skip == 1
            DIFF(1,ii) = nan; % record Diff under S_j(ii)
            RMSE(1,ii)= nan; % record rmse under S_j(ii)
            
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j_I(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
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
    disp(['     Parameterised OCV is not monotonic at SOC = ',num2str(S_j_I(indexOCV+1)), '    --> You MUST consider smoothing OCV.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%%% outlier detection for R0
Outlier_ind = isoutlier(R0,'movmedian',3);  % Hampel filter --> necessary but not sufficient
if any (Outlier_ind)
    Outlier_SOC = S_j_I(Outlier_ind);
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
    disp(['     Possible outlier(s) detected for R0 at SOC = ',num2str(S_j_I(Outlier_ind)), '    --> You MAY consider smoothing R0.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%% Final plots

cd Results
cd(fileNAME)
cd(newfile)

% Plots of OCV, R0, Ri, Ci, tau_i with SOC
f2=figure; plot (S_j_I, OCV, 'LineWidth',3); hold on; scatter(S_j_I, OCV,150,'k.'); xlabel('SOC'); ylabel('OCV (V)'); title (['Parameterised OCV with SOC @ ',num2str(Temps_depended(iiii)),'degC']); savefig(['Parameterised OCV @ ',num2str(Temps_depended(iiii)),'degC.fig']); 
f3=figure; plot (S_j_I, R0, 'LineWidth',3); hold on; scatter(S_j_I, R0,150,'k.'); xlabel('SOC'); ylabel('R0 (Ohm)'); title (['Parameterised R0 with SOC @ ',num2str(Temps_depended(iiii)),'degC']); savefig(['Parameterised R0 @ ',num2str(Temps_depended(iiii)),'degC.fig']); 
f4=figure; plot (S_j_I, R, 'LineWidth',3); hold on; scatter(S_j_I, R,150,'k.','HandleVisibility','off'); xlabel('SOC'); ylabel('Ri (Ohm)'); title (['Parameterised Ri with SOC @ ',num2str(Temps_depended(iiii)),'degC']); savefig(['Parameterised Ri @ ',num2str(Temps_depended(iiii)),'degC.fig']); 
f5=figure; plot (S_j_I, Ci, 'LineWidth',3); hold on; scatter(S_j_I, Ci,150,'k.','HandleVisibility','off'); xlabel('SOC'); ylabel('Ci (F)'); title (['Parameterised Ci with SOC @ ',num2str(Temps_depended(iiii)),'degC']); savefig(['Parameterised Ci @ ',num2str(Temps_depended(iiii)),'degC.fig']); 
f6=figure; plot (S_j_I, tau, 'LineWidth',3); hold on; scatter(S_j_I, tau,150,'k.','HandleVisibility','off'); xlabel('SOC'); ylabel('tau_i'); title (['Parameterised tau_i with SOC @ ',num2str(Temps_depended(iiii)),'degC']); savefig(['Parameterised tau_i @ ',num2str(Temps_depended(iiii)),'degC.fig']); 

cd ..
cd ..
cd ..

% Overall horizon - Experiment VS Paramterised Model
    
    if working_mode == 3 || working_mode == 4 % discharge  
      t_plot = t; I_plot = I; V_plot = V; T_plot = T; S_plot = S;
    else
      t_plot = t_charge; I_plot = I_charge; V_plot = V_charge; T_plot = T_charge; S_plot = S_charge;
    end
    
    plotit = [];
    count = 0;
      for kk = 1:numel(t_plot)
        indexTTT = ((Temps_depended(iiii) - Temps_sensitivity) <= T_plot{kk}) + (T_plot{kk} <= (Temps_depended(iiii) + Temps_sensitivity)) == 2; 
        
            if any(indexTTT)
                    count = count + 1;   
                    plotit.I {count,1} = I_plot{kk}(indexTTT);
                    plotit.t {count,1} = t_plot{kk}(indexTTT);
                    plotit.T {count,1} = T_plot{kk}(indexTTT);
                    plotit.V {count,1} = V_plot{kk}(indexTTT);
                    plotit.S {count,1} = S_plot{kk}(indexTTT); 
            end  
      end
    
    % filter out small piece of data  
    indxxx = cellfun(@length,plotit.T)>= (NRC * 3 + 2) * 2;
    plotit.I = plotit.I(indxxx);
    plotit.t = plotit.t(indxxx); 
    plotit.T = plotit.T(indxxx); 
    plotit.V = plotit.V(indxxx); 
    plotit.S = plotit.S(indxxx); 

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
   
    xlabel('Time (s)'); ylabel('Voltage (V)'); title (['Experiment VS Parameterised Model @ ',num2str(Temps_depended(iiii)),'degC - Global horizon']); 
        
      for jjjj = 1 : length(plotit.I)
        [Vhat, t_useful] = fullsolve(plotit.t{jjjj},plotit.I{jjjj},plotit.T{jjjj},plotit.S{jjjj},S_j_I,OCV,R0,R,tau,flags,NRC);
        if isnan(Vhat)
            DIFF_global{jjjj,1} = nan;
        else
            aaa= find (N(:,1) == t_useful(1)); bbb= find (N(:,1) == t_useful(end));
            if working_mode == 3 || working_mode == 1% anode
                if anode_voltage_sign < 0
                DIFF_global{jjjj,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
                else
                DIFF_global{jjjj,1} = Vhat + N(aaa:bbb,3);  % Global difference of segment kk
                end
            else
                DIFF_global{jjjj,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
            end  
            plot(t_useful,Vhat,'LineWidth',1); % Parameterised model data
        end
      end

      for i = 2: 1000
        leg_str{i} = ['Model result ',num2str(i-1)];  
      end
      leg_str{1}='Experiment'; 
      legend(leg_str); % legends for model results
    
    cd Results
    cd(fileNAME)
    cd(newfile)
    
    savefig(['Global comparison - Experiment VS Paramterised Model @ ',num2str(Temps_depended(iiii)),'degC.fig']);

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
LookUpTable_AllInOne = zeros(length(S_j_I),12)*nan;
LookUpTable_AllInOne (:,1) = S_j_I';
LookUpTable_AllInOne (:,2) = OCV';
LookUpTable_AllInOne (:,3) = R0';
[good,~] = size(R); [luck,~] = size(tau);
LookUpTable_AllInOne (:,4:3+good) = R';
LookUpTable_AllInOne (:,7:6+luck) = Ci';
LookUpTable_AllInOne (:,10:9+luck) = tau';

%%%%%%%%%%%% for overall collection and comparison %%%%%%%%%%%%%
TEMP_DADDY(iiii) = Temps_depended(iiii);
LUT_DADDY{iiii} = LookUpTable_AllInOne;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LookUpTable_AllInOne = mat2cell(LookUpTable_AllInOne,ones(1,length(S_j_I)),[1 1 1 1 1 1 1 1 1 1 1 1]);
title = {'SOC','OCV','R0','R1','R2','R3','C1','C2','C3','tau1','tau2','tau3'};
LookUpTable_AllInOne = [title;LookUpTable_AllInOne];

% writecell(LookUpTable_AllInOne,'LookUpTable_AllInOne.xlsx');

if flags.R0only.tf == 0
    disp(['All done @ ',num2str(Temps_depended(iiii)),' degC. The local and global root-mean-square errors (with respect to voltage) are ',num2str(rms(DIFF_vector)),' and ',num2str(RMSE_GLOBAL),'.'])
end
disp(['Now you can find parameterisation results in the folder \Results\',fileNAME,'\',newfile])
disp('***************************************************************************************************')
% figures poped up before




clearvars -except TEMP_DADDY LUT_DADDY plotit working_mode Ccap data DIFF DIFF_vector DIFF_global RMSE_GLOBAL flags I I_ I_charge N NRC OCV R Ci R0 RMSE S_real S S_ S_charge  S_interval S_j bounds t T t_ T_ t_charge T_charge tau tau_recommended V V_ V_charge datasheet_name fixedtau_realy_tf Nom_capacity OCVsheet_name OCVtable_tf RCpair_number SOC_window voltage_up_limit voltage_low_limit fileNAME Temps_depended Temps_sensitivity dataTTT S_j_I 
evalin('base', 'save(''INPUT.mat'')');
save OUTPUT

cd .. 
cd ..
cd ..

end

%% Post-process for all currents

PostProcess4Temp

%% Nothing can last forever
disp(['All done @ all temperatures. Now you can find parameterisation results in the folder \Results\',fileNAME])

diary off
cd ..
cd ..

toc
