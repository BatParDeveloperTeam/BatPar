function [OCV,R0,R,tau,DIFF_,skip_] = easyECMfit(t,I,V,T,S,NRC,S_j,SOCTab,OCVTab,SOCSOCTab,R0R0Tab,flags,Ccap,working_mode,R,tau) % [OCV,R0,R,tau_u,tau_l,initdata,Vhat]

    warning('off', 'all');
    
    skip_ = 0;

    N = numel(t);                         % number of segments in this window
    
if N == 0 % no data
    OCV = nan;
    R0 = nan;
    R = nan(NRC,1); 
    tau = nan(NRC,1);
    DIFF_ = nan;
    skip_ = 1;
    return
end

%% Calculate R0 from jumps by dV/dI
if  flags.taufix.tf  == 1  % Stage 2 enforcing; R0 and tau were already solved in Stage 1

    R0 = flags.R0solved.val;

else                               % Stage 1 enforcing;
    
    
    if (flags.SOCR0table.tf == 1)    % the SOC-R0 table exists
        if S_j >= min(SOCSOCTab) && S_j <= max(SOCSOCTab)
        R0  =  interp1(SOCSOCTab,R0R0Tab,S_j);
        end
    else
    
        R0_ini = cell(N,1); % initialise R0
    
        %   Tjumps = zeros(N,2); % please ignore this line - it is for T dependence
        
        for kk = 1:N  % for each segment%
            
        %     index0 = find (V{kk} > str2num(evalin('base','voltage_low_limit'))); % filter out noises that are below lower limit (voltage)
        %     I_R0 = I{kk}(index0);
        %     V_R0 = V{kk}(index0);
            
            I_R0 = I{kk};
            V_R0 = V{kk};
        
            dI = diff(I_R0);
            dV = diff(V_R0);
            
            if working_mode == 3 || working_mode == 4  % discharge
               if flags.pulsehead4R0.tf == 1
                % keep dI that are larger than 0.3C (This is assuming test noise is samller than 0.3C)
                    %index = find(abs(dI) > Ccap/3600 * 0.3); % both rising and dropping egdes considered 
                    index = find (dI >  Ccap/3600 * 0.28); % only rising egdes considered
               else 
                    index = find (dI < - Ccap/3600 * 0.28); % only dropping egdes considered
               end 
        
            else              % charge
               if flags.pulsehead4R0.tf == 1
                    % index = find(abs(dI) > Ccap/3600 * 0.3); % both rising and dropping egdes considered 
                    index = find (dI < - Ccap/3600 * 0.28); % only rising egdes considered
               else
                    index = find (dI > Ccap/3600 * 0.28); % only dropping egdes considered
               end
            end
        
            R0_ini{kk} = abs(dV(index)./dI(index));
        
        end 
            
            R0_ini = cell2mat(R0_ini);  % cell to matrix
            R0_ini = R0_ini(R0_ini ~= 0);    % filter out zero
            R0_ini = R0_ini(~isnan(R0_ini)); % filter out NaN
        
            R0 = mean(rmoutliers(R0_ini)); % filter out outliers and then average
            % R0 = mean(R0_ini); % filter out outliers and then average

    end
end 

if isnan(R0)
      OCV = nan;
      R0 = nan;
      R = nan(NRC,1); 
      tau = nan(NRC,1);
      DIFF_ = nan;
      skip_ = 1;
      return
end

if flags.R0only.tf == 1  % Parameterisation of R0 ONLY
      OCV = nan;
      R = nan(NRC,1); 
      tau = nan(NRC,1);
      DIFF_ = nan;
      return
end 
    
%% Formulate an optimisation problem to fit in fmincon - to solove OCV, Ri and Tau_i 


% OCV1    OCV2      Ri         tau_i            V0_i
% 1        2      3~2+NRC   3+NRC~2+2*NRC    3+2*NRC~2+2*NRC+N*NRC

% fmincon: https://uk.mathworks.com/help/optim/ug/fmincon.html
% [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

% param0 is 'x0'- intial values of unkonwns; 
param0 = zeros(2+(2+N)*NRC,1);  % 2{OCV} + 2{Ri & tau_i} *NRC + N{V0_i of each segment} *NRC
param0(3+NRC:2+2*NRC) = 5;      % tau_i=5
param0(3:2+NRC) = 0.001;      % R_i=0.001

% declaration of linear constraints {A,b,Aeq,beq} and options; no nonlinear constraints {nonlcon}
A = []; b = []; Aeq = []; beq = []; options = []; 

% lowerbound - lb
lb = -Inf*ones(2+(2+N)*NRC,1);  % initialisation of lowerbound; -Inf means no lowerbound
lb(3:2+NRC) = 0; % Ri > 0
lb(3+NRC:2+2*NRC) = 0; % tau_i > 0
% lb(3+2*NRC:end) = 0; % V0_i > 0   

% upperbound - ub
ub = Inf*ones(2+(2+N)*NRC,1);  % initialisation of upperbound; Inf means no lowerbound
ub(3:2+NRC) = 0.3; % Baseline : for most cells, its Ri <= 0.5
ub(3+NRC) = 1000; %tau_1 <= 100 %10
if NRC ==2
ub(3+NRC) = 50; %tau_1 <= 100 %50
ub(4+NRC) = 1500; % tau_2  <= 1500 
end
if NRC ==3
ub(3+NRC) = 10; %tau_1 <= 10 %10
ub(4+NRC) = 100; % tau_2  <= 100
ub(5+NRC) = 1000; % tau_3  <= 1000 
end

% if we relimit taui and Ri
if (flags.Ritaui_relimit.tf == 1)
    tau1_limit = str2num(evalin('base','tau1_limit'));
    tau2_limit = str2num(evalin('base','tau2_limit'));
    tau3_limit = str2num(evalin('base','tau3_limit'));
    Ri_limit = str2num(evalin('base','Ri_limit'));

    if isempty(tau1_limit) == 0
        ub(3+NRC) = tau1_limit;
    end 
    if isempty(tau2_limit) == 0
        ub(4+NRC) = tau2_limit;
    end 
    if isempty(tau3_limit) == 0
        ub(5+NRC) = tau3_limit;
    end 
    if isempty(Ri_limit) == 0
        ub(3:2+NRC) = Ri_limit;
    end 
end 

% if we are fixing tau (and also OCV) from stage 1
 if (flags.taufix.tf == 1)
 lb(3+NRC:2+2*NRC) = flags.taufix.val; % fix the lowerbound and upperbound the same at flags.taufix.val
 ub(3+NRC:2+2*NRC) = flags.taufix.val; 
 end

% the SOC-OCV table exists
if (flags.SOCOCVtable.tf == 1)  
    if S_j >= min(SOCTab) && S_j <= max(SOCTab)
    OCVbound  =  interp1(SOCTab,OCVTab,S_j);
    ub(1) = 0;
    lb(1) = 0;
    ub(2) = OCVbound;
    lb(2) = OCVbound;
    end
end

% (22 Apr 2022) fix the problem of inaccurate OCV param at 100% SOC discharge and at 0% SOC charge
if flags.SOCOCVtable.tf == 0  
    if ((working_mode == 3)||(working_mode == 4)) && (S_j == 1)  % dishcharge and S_j == 1
    ub(1) = 0;
    lb(1) = 0;
    ub(2) = max(cellfun(@max, V));
    lb(2) = max(cellfun(@max, V));
    % (22 Aug 2022) fix the problem of inaccurate R1 R2 param at 100% SOC discharge 
        for hhh = 1: NRC
            ub(3+hhh-1) = prctile (R(hhh,:),75);  % ub(3+hhh-1) = prctile (R(hhh,end-2:end),75);   
            ub(3+NRC+hhh-1) = prctile (tau(hhh,:),75);  % ub(3+NRC+hhh-1) = prctile (tau(hhh,end-2:end),75);  
        end
    end
    
    if ((working_mode == 1)||(working_mode == 2)) && (S_j == 0)  % charge and S_j == 0
    ub(1) = 0;
    lb(1) = 0;
    ub(2) = min(cellfun(@min, V));
    lb(2) = min(cellfun(@min, V));
    end
end

% Aeq * x = beq --> OCV-[V0_1 + V0_2+ V0_3] = V(1)+I(1)*R0   Use N constraints as initial conditions, important to convergence!!!
Aeq = zeros(N,2+(2+N)*NRC); % N by (1+(2+N)*NRC) matrix --> N segments(constraints), (1+(2+N)*NRC) rows in x
beq = zeros(N,1);           % N by 1 vector
for kk = 1:N   % for each segment(cosntraint)
    Aeq(kk,1) = S{kk}(1); % OCV1 * SOC
    Aeq(kk,2) = 1; % OCV2 * 1    
    Aeq(kk,3+(1+kk)*NRC:2+(2+kk)*NRC) = -1;  % [V0_1 + V0_2+ V0_3] * -1
    if (flags.R0temp.tf == 1)  % 
        beq(kk) = V{kk}(1) + I{kk}(1)*R0(T{kk}(1)); %
    else
        beq(kk) = V{kk}(1) + I{kk}(1)*R0;  % V(1)+I(1)*R0
    end
end

%  A * x <= b --> 
%  tau_1 - tau_2 <= 0, tau_2 - tau_3 <= 0, (NRC-1) constraints
%  -Vmax < V0_1 + V0_2+ V0_3 < Vmax,  N*2 contraints
%  Vmin < OCV1*S_j+OCV2 < Vmax,  2 constraints
    if (flags.SOCOCVtable.tf == 1)
        A = zeros(NRC-1+2*N,2+(2+N)*NRC); % NRC-1 by 2+(2+N)*NRC matrix --> NRC-1 constraints, (2+(2+N)*NRC) rows in x
        b = zeros(NRC-1+2*N,1);           % NRC-1 by 1 vector
    else
        A = zeros(NRC+2*N+1,2+(2+N)*NRC); % NRC+1 by 2+(2+N)*NRC matrix --> NRC+1 constraints, (2+(2+N)*NRC) rows in x
        b = zeros(NRC+2*N+1,1);           % NRC+1 by 1 vector
    end
    
    % constraints for tau
    for kk = 1:NRC-1              
        A(kk,2+NRC+kk) = 1;       
        A(kk,3+NRC+kk) = -1;      % if kk=1 --> tau_1 - tau_2 <= 0   if kk=2 --> tau_2 - tau_3 <= 0
    end
    
    % constraints for V0_i 
    for jj = 1:N
        for tao = 1:NRC
        A(NRC-1+jj,2+2*NRC+tao+(jj-1)*NRC) = 1;
        A(NRC-1+N+jj,2+2*NRC+tao+(jj-1)*NRC) = -1;
        end
        if flags.cathode_mode.tf == 0 && flags.anode_mode.tf == 0  % full cell
        b(NRC-1+jj) = str2num(evalin('base','voltage_up_limit')); % V0_1 + V0_2+ V0_3 <= 4.2
        b(NRC-1+N+jj) = str2num(evalin('base','voltage_up_limit')); % -[V0_1 + V0_2+ V0_3] <= 4.2 
        else
            if flags.cathode_mode.tf == 1 % cathode
            b(NRC-1+jj) = str2num(evalin('base','voltage_up_limit')) * 1.2; % a safety factor of 1.2 is placed for cathode
            b(NRC-1+N+jj) = str2num(evalin('base','voltage_up_limit')) * 1.2; % a safety factor of 1.2 is placed for cathode
            else  % anode
            b(NRC-1+jj) = str2num(evalin('base','voltage_up_limit')) * 0.2; % a safety factor of 0.2 is placed for cathode
            b(NRC-1+N+jj) = str2num(evalin('base','voltage_up_limit')) * 0.2; % a safety factor of 0.2 is placed for cathode    
            end 
        end
    end

    % constraints for OCV 
    if (flags.SOCOCVtable.tf == 0)
    A(NRC+2*N,1) = S_j;
    A(NRC+2*N,2) = 1;
    if flags.cathode_mode.tf == 0 && flags.anode_mode.tf == 0 % full cell
    b(NRC+2*N) = str2num(evalin('base','voltage_up_limit'));  % OCV1*S_j + OCV2 <= 4.2
    else
        if flags.cathode_mode.tf == 1 % cathode 
        b(NRC+2*N) = str2num(evalin('base','voltage_up_limit')) * 1.2; % a safety factor of 1.2 is placed for cathode  
        else % anode
        b(NRC+2*N) = 0; % Voc < 0  
        end
    end 
        if working_mode == 2  % full cell/cathode charge
            b(NRC+2*N) = min(cellfun(@max, V));    % OCV1*S_j + OCV2 <= min V(S_j(ii+1))   OCV(i)<OCV(ii+1)<V(ii+1)
        end

    A(NRC+2*N+1,1) = -S_j;
    A(NRC+2*N+1,2) = -1;
    b(NRC+2*N+1) = - str2num(evalin('base','voltage_low_limit'));  % -OCV1*S_j - OCV2 <= - 2.7
        if working_mode == 4  % full cell/cathode discharge
            b(NRC+2*N+1) = - max(cellfun(@min, V));  % -OCV1*S_j - OCV2 <= - max V(S_j(ii-1))  OCV(i)>OCV(ii-1)>V(ii-1)
        end 
        if working_mode == 1  % anode - charge
            b(NRC+2*N+1) = str2num(evalin('base','voltage_up_limit')) * 0.2 ; % -OCV1*S_j - OCV2 <= 0.2* Vmax  (OCV1*S_j + OCV2 >= -0.2*Vmax)
        end
        if working_mode == 3  % anode - discharge
            % b(NRC+2*N+1) = str2num(evalin('base','voltage_up_limit')) * 0.2 ; % -OCV1*S_j - OCV2 <= 0.2* Vmax  (OCV1*S_j + OCV2 >= -0.2*Vmax)
            b(NRC+2*N+1) = - max(cellfun(@min, V)); % -OCV1*S_j - OCV2 <= -V   V(ii-1) < OCV(ii-1) < OCV(ii)  -OCV(ii) < -OCV(ii-1) < -V(ii-1)
            
        end
    end 

%% Sovle formulated optmisation problem by fmincon

% options = optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',NRC*1500,'MaxIterations',NRC*500,'OptimalityTolerance',1e-20,'StepTolerance', 1e-20, 'FunctionTolerance', 1e-20);
options = optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',NRC*1500,'MaxIterations',NRC*500,'OptimalityTolerance',1e-20,'StepTolerance', 1e-20, 'FunctionTolerance', 1e-20);

% default values: MaxFunctionEvaluations:3000; MaxIterations:1000; OptimalityTolerance:1e-6; StepTolerance:1e-10; FunctionTolerance:1e-6
% In general, you can improve accuracy of params by increasing -ations or decreasing tolerances, but it is not always worth doing - tiny improvements VS significant computational loads, and sometimes overfitting!!!

x_solved = fmincon(@(param) F_easyECMfit(t,I,V,T,S,NRC,R0,flags,param),param0,A,b,Aeq,beq,lb,ub,[],options); % x_solved - column vector

% split up params to output    
OCV1 = x_solved(1); OCV2 = x_solved(2); OCV = OCV1*S_j + OCV2; % optimised value of the OCV
R = x_solved(3:2+NRC); tau = x_solved(3+NRC:2+2*NRC);          % optimised value of Ri and tau_i
V0_i = x_solved(3+2*NRC:end); initdata = transpose(reshape(V0_i,[NRC N]));

%% difference of voltage between experiment and parameterised model

DIFF_ = [];  % Intialise a vector connector 
for kk = 1:N
        Diff_Exp_Mod = localsolve(t{kk},I{kk},T{kk},S{kk},OCV1,OCV2,R0,R,tau,initdata(kk,:),flags)-V{kk}; % difference between experimental and model voltage
        Diff_Exp_Mod = Diff_Exp_Mod(~isnan(Diff_Exp_Mod)); % filter out NaN       
        DIFF_ = cat(1,DIFF_,Diff_Exp_Mod); % Connect two column vectors; '1' menas along column direction 
end 
%% Plots of voltage of experiment and parameterised model

if (flags.allplots.tf == 1)
    figure; 
    xlabel('Relative time (s)');
    ylabel('Voltage (V)');
    title (['Experiment VS Paramterised Model - Voltage at around SOC = ', num2str(S_j),]);
    hold on;
    for kk = 1:N
        plot(t{kk}-t{kk}(1),V{kk}); legend (['Experiment data ',num2str(kk),]);
        plot(t{kk}-t{kk}(1),localsolve(t{kk},I{kk},T{kk},S{kk},OCV1,OCV2,R0,R,tau,initdata(kk,:),flags)); legend (['Model data ',num2str(kk),]);
    end
    
    for iii = 1:N
    leg_string{2*iii-1} = ['Experiment data ',num2str(iii)];  
    leg_string{2*iii} = ['Model data ',num2str(iii)];
    end
    legend(leg_string); % legends for experiment and model data
       
    drawnow;  % normally plot commands incorporated in fuctions would only be executed after running over the main programme.drawnow is here to enfore plot commands immediately
end
end
