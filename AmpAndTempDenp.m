

[dataIT, SOC_steps] = data_IandTdependent(data, Currents_depended, Currents_sensitivity, Temps_depended, Temps_sensitivity, S_j , NRC);

if flags.R0only.tf == 0
    if flags.taufixed_Stage2.tf == 0
    disp(['Initialisation done - Paramterisation (with VARIABLE tau and I dependence and T dependence) in progress, including ',num2str(length(Currents_depended)*length(Temps_depended)), ' PHASES (', num2str(length(Currents_depended)),' CURRENTS by ',num2str(length(Temps_depended)),' TEMPERATURES) --> Each PHASE further includes ~',num2str(max(max(SOC_steps))),' STEPS.']) 
    else
    disp(['Initialisation done - Paramterisation (with CONSTANT tau and I dependence and T dependence) in progress, including ',num2str(length(Currents_depended)*length(Temps_depended)), ' PHASES (', num2str(length(Currents_depended)),' CURRENTS by ',num2str(length(Temps_depended)),' TEMPERATURES) --> Each PHASE further includes 2 STAGES --> Each STAGE further includes ~',num2str(max(max(SOC_steps))),' STEPS.']) 
    end
    disp('***************************************************************************************************')
else 
    disp(['Initialisation done - Paramterisation (with I dependence and T dependence) in progress, including ',num2str(length(Currents_depended)*length(Temps_depended)), ' PHASES (', num2str(length(Currents_depended)),' CURRENTS by ',num2str(length(Temps_depended)),' TEMPERATURES) --> Each PHASE further includes ~',num2str(max(max(SOC_steps))),' STEPS.']) 
    disp('***************************************************************************************************')
end 

for iiii = 1:length(Currents_depended)
    for jjjj = 1:length(Temps_depended)
C_rate_round = roundn(Currents_depended(iiii)/(Ccap/3600),-2);  % rounded C rate
phase_number = (iiii - 1)*length(Temps_depended) + jjjj;        % Phase number

disp(['-> Parameterisation taking place at PHASE ',num2str(phase_number),' of ',num2str(length(Currents_depended)*length(Temps_depended)),' -- ', num2str(C_rate_round),'C (',num2str(Currents_depended(iiii)),'Amp) & ',num2str(Temps_depended(jjjj)),'degC']) 

if working_mode == 3 || working_mode == 4
    f111=figure;
    hold on
    for kk = 1:numel(I) % discharge segment
        indexII = currentIndexing (I{kk}, Currents_depended(iiii), Currents_sensitivity);
        indexTT = ((Temps_depended(jjjj) - Temps_sensitivity) <= T{kk}) + (T{kk} <= (Temps_depended(jjjj) + Temps_sensitivity)) == 2; 
        indexIITT = indexII + indexTT == 2;        
        plot(S{kk}(indexIITT),V{kk}(indexIITT))
        xlabel('SOC');
        ylabel('Voltage (V)');
        title (['Data extraction (discharge) - Voltage VS SOC @ ',num2str(C_rate_round),'C & ', num2str(Temps_depended(jjjj)),'degC']);
    end
    
cd Results
cd(fileNAME)
newfile = [fileNAME,'_at_',num2str(C_rate_round),'C_and_',num2str(Temps_depended(jjjj)),'degC'];
mkdir(newfile)
cd(newfile)
    
savefig(['Data extraction of discharge segments @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig'])

else
    f222=figure;
    hold on
    for kk = 1:numel(I_charge) % discharge segment
        indexII = currentIndexing (I_charge{kk}, Currents_depended(iiii), Currents_sensitivity);
        indexTT = ((Temps_depended(jjjj) - Temps_sensitivity) <= T_charge{kk}) + (T_charge{kk} <= (Temps_depended(jjjj) + Temps_sensitivity)) == 2; 
        indexIITT = indexII + indexTT == 2;    
        plot(S_charge{kk}(indexIITT),V_charge{kk}(indexIITT))
        xlabel('SOC');
        ylabel('Voltage (V)');
        title (['Data extraction (charge) - Voltage VS SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']);
    end

cd Results
cd(fileNAME)
newfile = [fileNAME,'_at_',num2str(C_rate_round),'C_and_',num2str(Temps_depended(jjjj)),'degC'];
mkdir(newfile)
cd(newfile)    
savefig(['Data extraction of charge segments @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig'])
end


cd ..
cd ..
cd ..


if flags.taufixed_Stage2.tf == 1 && flags.R0only.tf == 0
disp('--> Proceeding to Stage 1 - solving parameters with VARIABLE tau.') 
end

S_j_IT = dataIT.S_j_IT{iiii, jjjj};  % the real SOC at the current current and temp

% Initialisation of unkonws (variables to be parameterised)
OCV = zeros(1,numel(S_j_IT));   
R0  = zeros(1+flags.R0temp.tf,numel(S_j_IT));
R   = zeros(NRC,numel(S_j_IT));
tau = zeros(NRC,numel(S_j_IT));

%---------- Initialisation for the calculation of RMSE
L_temp = [];   % wipe out values assigned in previous operation
for tao = 1:numel(S_j_IT)
L_temp(tao) = sum(cellfun(@length, dataIT.t{iiii,jjjj}{1,tao})); % sum length of all martrixes of one cell of data.t
end
L_dim = max(L_temp); % figure out the max(sum length of all martrixes of each cell of data.t)

DIFF = zeros(L_dim,numel(S_j_IT))*nan; % the matrix records volatge difference between experiment and parameterised model

RMSE = zeros(1,numel(S_j_IT)); % the vector records the root-mean-square error between experiment and model at different S-j
%----------

%%% The show starts
flags.taufix.tf  = 0;  % provisionally disabling fixed tau functions
for ii = 1:numel(S_j_IT)
    
    disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j_IT)),', parameterising at SOC ', num2str(S_j_IT(ii)),])  % display the progress of parameterisation
                       
    t_ = dataIT.t{iiii,jjjj}{ii};  
    I_ = dataIT.I{iiii,jjjj}{ii}; 
    V_ = dataIT.V{iiii,jjjj}{ii};  % pass data
    T_ = dataIT.T{iiii,jjjj}{ii};
    S_ = dataIT.S{iiii,jjjj}{ii};
    
    [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j_IT(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
    % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
    
    if Skip == 1
        DIFF(1,ii) = nan; % record Diff under S_j(ii)
        RMSE(1,ii)= nan; % record rmse under S_j(ii)
        
        disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j_IT(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
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
    
    for ii = 1:numel(S_j_IT)
        
        disp(['     Step ',num2str(ii),' of ',num2str(numel(S_j_IT)),', parameterising at SOC ', num2str(S_j_IT(ii)),])  % display the progress of parameterisation
                           
         t_ = dataIT.t{iiii,jjjj}{ii};  
         I_ = dataIT.I{iiii,jjjj}{ii}; 
         V_ = dataIT.V{iiii,jjjj}{ii};  % pass data
         T_ = dataIT.T{iiii,jjjj}{ii};
         S_ = dataIT.S{iiii,jjjj}{ii};
         flags.R0solved.val  = R0(ii); % Solved R0 from Stage 1

        [OCV(ii),R0(ii),R(:,ii),tau(:,ii), Diff, Skip ] = easyECMfit(t_,I_,V_,T_,S_,NRC,S_j_IT(ii),SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau); % the expression to optimise OCV, R0, R_i, tau_i
        % OCV - 1 by numel(S_j), R - NRC by numel(S_j), tau - NRC by numel(S_j)
        if Skip == 1
            DIFF(1,ii) = nan; % record Diff under S_j(ii)
            RMSE(1,ii)= nan; % record rmse under S_j(ii)
            
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp(['     Step ',num2str(ii),' failed. --> Parameterisation cannot happen at SOC = ',num2str(S_j_IT(ii)),' because of insufficient experimental data. The codes have to jump over this SOC. If you keep seeing this message, please make the SOC window size larger.'])
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
    disp(['     Parameterised OCV is not monotonic at SOC = ',num2str(S_j_IT(indexOCV+1)), '    --> You MUST consider smoothing OCV.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%%% outlier detection for R0
Outlier_ind = isoutlier(R0,'movmedian',3);  % Hampel filter --> necessary but not sufficient
if any (Outlier_ind)
    Outlier_SOC = S_j_IT(Outlier_ind);
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
    disp(['     Possible outlier(s) detected for R0 at SOC = ',num2str(S_j_IT(Outlier_ind)), '    --> You MAY consider smoothing R0.'])
    disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
end

%% Final plots

cd Results
cd(fileNAME)
cd(newfile)

% Plots of OCV, R0, Ri, Ci, tau_i with SOC
f2=figure; plot (S_j_IT, OCV, 'LineWidth',3); hold on; scatter(S_j_IT, OCV,150,'k.'); xlabel('SOC'); ylabel('OCV (V)'); title (['Parameterised OCV with SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']); savefig(['Parameterised OCV @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']); 
f3=figure; plot (S_j_IT, R0, 'LineWidth',3); hold on; scatter(S_j_IT, R0,150,'k.'); xlabel('SOC'); ylabel('R0 (Ohm)'); title (['Parameterised R0 with SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']); savefig(['Parameterised R0 @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']); 
f4=figure; plot (S_j_IT, R, 'LineWidth',3); xlabel('SOC'); ylabel('Ri (Ohm)'); title (['Parameterised Ri with SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']); hold on; scatter(S_j_IT, R,150,'k.','HandleVisibility','off'); savefig(['Parameterised Ri @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']); 
f5=figure; plot (S_j_IT, Ci, 'LineWidth',3); xlabel('SOC'); ylabel('Ci (F)'); title (['Parameterised Ci with SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']); hold on; scatter(S_j_IT, Ci,150,'k.','HandleVisibility','off'); savefig(['Parameterised Ci @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']); 
f6=figure; plot (S_j_IT, tau, 'LineWidth',3); xlabel('SOC'); ylabel('tau_i'); title (['Parameterised tau_i with SOC @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC']); hold on; scatter(S_j_IT, tau,150,'k.','HandleVisibility','off'); savefig(['Parameterised tau_i @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']); 

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
        
          indexTTT = ((Temps_depended(jjjj) - Temps_sensitivity) <= T_plot{kk}) + (T_plot{kk} <= (Temps_depended(jjjj) + Temps_sensitivity)) == 2;
          T_plot{kk} = T_plot{kk}(indexTTT);
          I_plot{kk} = I_plot{kk}(indexTTT);
          V_plot{kk} = V_plot{kk}(indexTTT);
          S_plot{kk} = S_plot{kk}(indexTTT);
          t_plot{kk} = t_plot{kk}(indexTTT);
 
          if ~isempty(I_plot{kk})
          indexIII = currentIndexing (I_plot{kk}, Currents_depended(iiii), Currents_sensitivity);
        
            if any(indexIII)
                [~, minlocs] = findpeaks(-indexIII);
                if ~isempty(minlocs)
                    count = count + 1;  
                    indexIII_temp = indexIII;
                    indexIII_temp(minlocs(1) : end) = false;
                    plotit.I {count,1} = I_plot{kk}(indexIII_temp);
                    plotit.t {count,1} = t_plot{kk}(indexIII_temp);
                    plotit.T {count,1} = T_plot{kk}(indexIII_temp);
                    plotit.V {count,1} = V_plot{kk}(indexIII_temp);
                    plotit.S {count,1} = S_plot{kk}(indexIII_temp);
                    
                    count = count + 1;  
                    indexIII_temp = indexIII;
                    indexIII_temp(1 : minlocs(end)) = false;
                    plotit.I {count,1} = I_plot{kk}(indexIII_temp);
                    plotit.t {count,1} = t_plot{kk}(indexIII_temp);
                    plotit.T {count,1} = T_plot{kk}(indexIII_temp);
                    plotit.V {count,1} = V_plot{kk}(indexIII_temp);
                    plotit.S {count,1} = S_plot{kk}(indexIII_temp);
                    
                    if  length(minlocs) > 1
                        for hhh = 2: length(minlocs)
                            count = count + 1;
                            indexIII_temp = indexIII;
                            indexIII_temp(1 : minlocs(hhh - 1)) = false;
                            indexIII_temp(minlocs(hhh) : end) = false;
                            plotit.I {count,1} = I_plot{kk}(indexIII_temp);
                            plotit.t {count,1} = t_plot{kk}(indexIII_temp);
                            plotit.T {count,1} = T_plot{kk}(indexIII_temp);
                            plotit.V {count,1} = V_plot{kk}(indexIII_temp);
                            plotit.S {count,1} = S_plot{kk}(indexIII_temp);
                        end 
                    end 

                else
                    count = count + 1;   
                    plotit.I {count,1} = I_plot{kk}(indexIII);
                    plotit.t {count,1} = t_plot{kk}(indexIII);
                    plotit.T {count,1} = T_plot{kk}(indexIII);
                    plotit.V {count,1} = V_plot{kk}(indexIII);
                    plotit.S {count,1} = S_plot{kk}(indexIII);
                end 
            end
          end  
      end
    
    % filter out small piece of data  
    indxxx = cellfun(@length,plotit.I)>= (NRC * 3 + 2) * 2;
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
    
    xlabel('Time (s)'); ylabel('Voltage (V)'); title (['Experiment VS Parameterised Model @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC - Global horizon']); 
        
      for kkkk = 1 : length(plotit.I)
        [Vhat, t_useful] = fullsolve(plotit.t{kkkk},plotit.I{kkkk},plotit.T{kkkk},plotit.S{kkkk},S_j_IT,OCV,R0,R,tau,flags,NRC);
        if isnan(Vhat)
            DIFF_global{kkkk,1} = nan;
        else
            aaa= find (N(:,1) == t_useful(1)); bbb= find (N(:,1) == t_useful(end));
            if working_mode == 3 || working_mode == 1% anode
                if anode_voltage_sign < 0
                DIFF_global{kkkk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
                else
                DIFF_global{kkkk,1} = Vhat + N(aaa:bbb,3);  % Global difference of segment kk
                end
            else
                DIFF_global{kkkk,1} = Vhat - N(aaa:bbb,3);  % Global difference of segment kk
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
    
    savefig(['Global comparison - Experiment VS Paramterised Model @ ',num2str(C_rate_round),'C & ',num2str(Temps_depended(jjjj)),'degC.fig']);

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
LookUpTable_AllInOne = zeros(length(S_j_IT),12)*nan;
LookUpTable_AllInOne (:,1) = S_j_IT';
LookUpTable_AllInOne (:,2) = OCV';
LookUpTable_AllInOne (:,3) = R0';
[good,~] = size(R); [luck,~] = size(tau);
LookUpTable_AllInOne (:,4:3+good) = R';
LookUpTable_AllInOne (:,7:6+luck) = Ci';
LookUpTable_AllInOne (:,10:9+luck) = tau';

%%%%%%%%%%%% for overall collection and comparison %%%%%%%%%%%%%
C_RATE_DADDY(iiii) = C_rate_round;
TEMP_DADDY(jjjj) = Temps_depended(jjjj);
LUT_DADDY{iiii,jjjj} = LookUpTable_AllInOne;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LookUpTable_AllInOne = mat2cell(LookUpTable_AllInOne,ones(1,length(S_j_IT)),[1 1 1 1 1 1 1 1 1 1 1 1]);
title = {'SOC','OCV','R0','R1','R2','R3','C1','C2','C3','tau1','tau2','tau3'};
LookUpTable_AllInOne = [title;LookUpTable_AllInOne];
writecell(LookUpTable_AllInOne,'LookUpTable_AllInOne.xlsx');

if flags.R0only.tf == 0
    disp(['All done @ ',num2str(C_rate_round),'C & ', num2str(Temps_depended(jjjj)),'degC. The local and global root-mean-square errors (with respect to voltage) are ',num2str(rms(DIFF_vector)),' and ',num2str(RMSE_GLOBAL),'.'])
end
disp(['Now you can find parameterisation results in the folder \Results\',fileNAME,'\',newfile])
disp('***************************************************************************************************')
% figures poped up before

clearvars -except iiii TEMP_DADDY C_RATE_DADDY LUT_DADDY phase_number plotit working_mode Ccap data DIFF DIFF_vector DIFF_global RMSE_GLOBAL flags I I_ I_charge N NRC OCV R Ci R0 RMSE S_real S S_ S_charge  S_interval S_j bounds t T t_ T_ t_charge T_charge tau tau_recommended V V_ V_charge datasheet_name fixedtau_realy_tf Nom_capacity OCVsheet_name OCVtable_tf RCpair_number SOC_window voltage_up_limit voltage_low_limit fileNAME Currents_depended Currents_sensitivity Temps_depended Temps_sensitivity dataIT S_j_IT C_rate_round 
evalin('base', 'save(''INPUT.mat'')');
save OUTPUT

cd .. 
cd ..
cd ..

    end
end

%% Post-process for all currents and temperatures

PostProcess4CrateandTemps

%% Nothing can last forever - except thoughtful codes
disp(['All done @ all C-rates and Temperatures. Now you can find parameterisation results in the folder \Results\',fileNAME])

diary off
cd ..
cd ..

toc
