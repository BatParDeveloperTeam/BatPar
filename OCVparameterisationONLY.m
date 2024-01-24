SOC_param_lowlimit = str2num(evalin('base','SOC_param_lowlimit'));
SOC_param_uplimit = str2num(evalin('base','SOC_param_uplimit'));

%% True OCV
if strcmp(flags.OCVtrue.tf,'yes')
    
    % Cropping SOC range   
    if isempty(SOC_param_lowlimit) == 0
        index_SOC = S_ >= SOC_param_lowlimit;
        S_ = S_(index_SOC);
        I_ = I_(index_SOC);
        V_ = V_(index_SOC);
        t_ = t_(index_SOC);
    end
    
    if isempty(SOC_param_uplimit) == 0
        index_SOC = S_ <= SOC_param_uplimit;
        S_ = S_(index_SOC);
        I_ = I_(index_SOC);
        V_ = V_(index_SOC);
        t_ = t_(index_SOC);
    end   
    
    disp(['SOC range considered in parameterisation is between ',num2str(min(S_)),' and ', num2str(max(S_)),'.'])
    disp('Initialisation done - true OCV paramterisation in progress')


    % OCV extraction
    j = 1; 
    for i = 1:length(I_)-1
         if (I_(i) == 0) && (I_(i+1) ~= 0)
            SOC_OCV (j,1) = S_ (i);
            SOC_OCV (j,2) = V_ (i);
            j = j + 1;
         end 
    end
    
    if I_(end) == 0
        SOC_OCV (j,1) = S_ (end);
        SOC_OCV (j,2) = V_ (end);
    end
    SOC_OCV = sortrows(SOC_OCV,1);

    % post-processing
    cd Results
    cd(fileNAME)
    
    f217=figure;
    plot (SOC_OCV (:,1), SOC_OCV (:,2)); hold on; scatter(SOC_OCV (:,1), SOC_OCV (:,2),150,'k.'); title ('Parameterised OCV with SOC'); savefig('Parameterised OCV.fig'); 
    xlabel('SOC'); ylabel('OCV (V)'); 
    
    LookUpTable_AllInOne = zeros(length(SOC_OCV (:,1)),2)*nan;
    LookUpTable_AllInOne (:,1) = SOC_OCV (:,1);
    LookUpTable_AllInOne (:,2) = SOC_OCV (:,2);
    LookUpTable_AllInOne = mat2cell(LookUpTable_AllInOne,ones(1,length(SOC_OCV (:,1))),[1 1]);
    title = {'SOC','OCV'};
    LookUpTable_AllInOne = [title;LookUpTable_AllInOne];
    writecell(LookUpTable_AllInOne,'LookUpTable_AllInOne.xlsx');
    disp(['Now you can find parameterisation results in the folder \Results\',fileNAME])
    
    evalin('base', 'save(''INPUT.mat'')');
    save OUTPUT

    toc
    diary off

    cd ..
    cd ..
end 


%% Pseudo OCV
if strcmp(flags.OCVpseudo.tf,'yes')
    
    % charge or discharge
    if strcmp(flags.chagre_mode.tf,'yes')
        S_ = S_charge;
        I_ = I_charge;
        V_ = V_charge;
        t_ = t_charge;
    else 
        S_ = S;
        I_ = I;
        V_ = V;
        t_ = t;
    end
    
    % Cropping SOC range
    if isempty(SOC_param_lowlimit) == 0
        for k = 1: numel(S_)
            index_SOC = S_{k,1} >= SOC_param_lowlimit;
            S_{k,1} = S_{k,1}(index_SOC);
            I_{k,1} = I_{k,1}(index_SOC);
            V_{k,1} = V_{k,1}(index_SOC);
            t_{k,1} = t_{k,1}(index_SOC);
        end 
    end

    if isempty(SOC_param_uplimit) == 0
        for k = 1: numel(S_)
            index_SOC = S_{k,1} <= SOC_param_uplimit;
            S_{k,1} = S_{k,1}(index_SOC);
            I_{k,1} = I_{k,1}(index_SOC);
            V_{k,1} = V_{k,1}(index_SOC);
            t_{k,1} = t_{k,1}(index_SOC);
        end 
    end
    
    % ascending sortation
    SOC_old = cell2mat(S_);
    OCV_old = cell2mat(V_);  
    [~,ind_,~] = unique(SOC_old); % (~:neglected  ind: index of unique rows)
    SOC_OCV_old(:,1) = SOC_old(ind_);
    SOC_OCV_old(:,2) = OCV_old(ind_);
    SOC_OCV_old = sortrows(SOC_OCV_old,1);
        
    disp(['SOC range considered in parameterisation is between ',num2str(min(SOC_OCV_old(:,1))),' and ', num2str(max(SOC_OCV_old(:,1))),'.'])
    disp('Initialisation done - pseudo OCV paramterisation in progress')

    % OCV interpolation
    SOC_new = min(SOC_OCV_old(:,1)) : 1e-4 : max(SOC_OCV_old(:,1));  % new SOC trace
    SOC_new = SOC_new';
    OCV_new = interp1(SOC_OCV_old(:,1), SOC_OCV_old(:,2), SOC_new);  % generate new OCV trace;


    % post-processing
    cd Results
    cd(fileNAME)
    
    f218=figure;
    plot (SOC_new, OCV_new); title ('Parameterised OCV with SOC'); savefig('Parameterised OCV.fig'); 
    xlabel('SOC'); ylabel('OCV (V)'); 
    
    LookUpTable_AllInOne = zeros(length(SOC_new),2)*nan;
    LookUpTable_AllInOne (:,1) = SOC_new;
    LookUpTable_AllInOne (:,2) = OCV_new;
    LookUpTable_AllInOne = mat2cell(LookUpTable_AllInOne,ones(1,length(SOC_new)),[1 1]);
    title = {'SOC','OCV'};
    LookUpTable_AllInOne = [title;LookUpTable_AllInOne];
    writecell(LookUpTable_AllInOne,'LookUpTable_AllInOne.xlsx');
    disp(['Now you can find parameterisation results in the folder \Results\',fileNAME])
    
    evalin('base', 'save(''INPUT.mat'')');
    save OUTPUT

    toc
    diary off

    cd ..
    cd ..

end