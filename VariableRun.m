clc

todolist = readcell('list.xlsx');  % read list.xlsx
[rows_, columns_] = size(todolist);

for iiiiii = 3: columns_
    if ismissing(todolist{5, iiiiii})
        column_dataend = iiiiii-1;    % the last datasheet specified
        break
    else
        column_dataend = columns_;
    end
end

%% Check for missing info
countMiss = 0;
MissingInfo = [];
for jjjjjj = 2:19 
    for kkkkkk = 3: column_dataend
        if ismissing(todolist{jjjjjj, kkkkkk})
            countMiss = countMiss + 1;
            MissingInfo (countMiss) = kkkkkk - 2;            
        end
    end
end
MissingInfo = unique(MissingInfo);

% if ~isempty (MissingInfo)
%     disp(['There are missing inputs in list.xlsx, specifically in Data ',num2str(MissingInfo),'. Please go double check !'])
%     return
% end 


disp(['Batch run started, using inputs defined in {list.xlsx} - ',num2str(column_dataend - 2),' datasheets found. Good luck!'])


%% loop - read inputs and run

assignin('base','todolist',todolist);
assignin('base','column_dataend',column_dataend);

for tttttttt = 3: column_dataend 
    clearvars -except tttttttt Files todolist column_dataend
    todolist = evalin('base','todolist');
    column_dataend = evalin('base','column_dataend');

    % full or half cell
    if  strcmp(todolist{2,tttttttt},'Full cell')    % full cell
        assignin('base','Mode_anode_tf','no');
        assignin('base','Mode_cathode_tf','no');
    else 
        if strcmp(todolist{2,tttttttt},'Cathode')     % Cathode
            assignin('base','Mode_anode_tf','no');
            assignin('base','Mode_cathode_tf','yes');           
        else                                           % Anode
            assignin('base','Mode_anode_tf','yes');   
            assignin('base','Mode_cathode_tf','no');
        end
    end
    
    % discharge or charge
    if  strcmp(todolist{3,tttttttt},'Charge')    % discharge
        assignin('base','Mode_charge_tf','yes');
    else 
        assignin('base','Mode_charge_tf','no');
    end    

    assignin('base','exp_path',todolist{4,tttttttt});
    
    assignin('base','datasheet_name',todolist{5,tttttttt});

    if ismissing(todolist{6, tttttttt})
        assignin('base','record_column',num2str([]));
    else
        assignin('base','record_column',num2str(todolist{6,tttttttt}));
    end 

    assignin('base','time_column',num2str(todolist{7,tttttttt}));
    
    assignin('base','current_column',num2str(todolist{8,tttttttt}));
    
    assignin('base','voltage_column',num2str(todolist{9,tttttttt}));
    
    assignin('base','temp_column',num2str(todolist{10,tttttttt}));

    % current sign
    if strcmp(todolist{11,tttttttt},'Below zero')
        assignin('base','Current_belwozero_tf','yes');
    else
        assignin('base','Current_belwozero_tf','no'); 
    end  
    
    assignin('base','Nom_capacity',num2str(todolist{12,tttttttt}));
    
    assignin('base','voltage_up_limit',num2str(todolist{13,tttttttt}));

    assignin('base','voltage_low_limit',num2str(todolist{14,tttttttt}));

    % SOC baseline is 1 or 0
    if strcmp(todolist{15,tttttttt},'Yes')
        assignin('base','MaxSOCis1','yes');
    else
        assignin('base','MaxSOCis1','no'); 
    end

    assignin('base','RCpair_number',num2str(todolist{16,tttttttt}));

    assignin('base','SOC_window',num2str(todolist{17,tttttttt}));

    % tau to be constant or variable
    if strcmp(todolist{18,tttttttt},'Constant')
        assignin('base','fixedtau_realy_tf','yes');
    else
        assignin('base','fixedtau_realy_tf','no'); 
    end
    
    % pulse head or end
    if strcmp(todolist{19,tttttttt},'Pulse head')
        assignin('base','Pulse_head_tf','yes');
    else
        assignin('base','Pulse_head_tf','no'); 
    end

    % SOC high limit
    if ismissing(todolist{21, tttttttt})
        assignin('base','SOC_param_lowlimit',num2str([]));
    else
        assignin('base','SOC_param_lowlimit',num2str(todolist{21,tttttttt}));
    end 
    
    % SOC low limit
    if ismissing(todolist{22, tttttttt})
        assignin('base','SOC_param_uplimit',num2str([]));
    else
        assignin('base','SOC_param_uplimit',num2str(todolist{22,tttttttt}));
    end 

    % datasheet cut
    if ismissing(todolist{23, tttttttt}) && ismissing(todolist{24, tttttttt})
        assignin('base','datasheet_cleanup_tf','no');
    else 
        assignin('base','datasheet_cleanup_tf','yes');
    end   
    if ismissing(todolist{23, tttttttt})
        assignin('base','row_start',num2str([]));
    else
        assignin('base','row_start',num2str(todolist{23,tttttttt}));
    end 
    if ismissing(todolist{24, tttttttt})
        assignin('base','row_end',num2str([]));
    else
        assignin('base','row_end',num2str(todolist{24,tttttttt}));
    end 

    % SOC-OCV table
    if ismissing(todolist{25, tttttttt})
        assignin('base','OCVtable_tf','no');
        assignin('base','OCVsheet_name',num2str([]));
    else
        assignin('base','OCVtable_tf','yes');
        assignin('base','OCVsheet_name',todolist{25, tttttttt});
    end 
    
    % SOC-R0 table
    if ismissing(todolist{26, tttttttt})
        assignin('base','R0table_tf','no');
        assignin('base','R0sheet_name',num2str([]));
    else
        assignin('base','R0table_tf','yes');
        assignin('base','R0sheet_name',todolist{26, tttttttt});
    end 

    % Parameterisation of OCV or R0 only
    if  strcmp(todolist{27,tttttttt},'None')    % No
        assignin('base','OCVpseudo_tf','no');
        assignin('base','OCVtrue_tf','no');   
        assignin('base','R0only_tf','no'); 
    else 
        if strcmp(todolist{27,tttttttt},'Pseudo OCV from constant current')     % Cathode
            assignin('base','OCVpseudo_tf','yes');
            assignin('base','OCVtrue_tf','no');   
            assignin('base','R0only_tf','no');         
        else                                           
            if strcmp(todolist{27,tttttttt},'True OCV from relaxations')     % Cathode
                assignin('base','OCVpseudo_tf','no');
                assignin('base','OCVtrue_tf','yes');   
                assignin('base','R0only_tf','no'); 
            else            
                assignin('base','OCVpseudo_tf','no');
                assignin('base','OCVtrue_tf','no');   
                assignin('base','R0only_tf','yes'); 
            end
        end
    end

    % Current dependence
    if ismissing(todolist{28, tttttttt}) || ismissing(todolist{29, tttttttt})
        assignin('base','CurrentDenp_tf','no');
    else 
        assignin('base','CurrentDenp_tf','yes');
    end   
    if ismissing(todolist{28, tttttttt})
        assignin('base','Currents_4_Denp',num2str([]));
    else
        assignin('base','Currents_4_Denp',num2str(todolist{28,tttttttt}));
    end 
    if ismissing(todolist{29, tttttttt})
        assignin('base','CurrentSensi_4_Denp',num2str([]));
    else
        assignin('base','CurrentSensi_4_Denp',num2str(todolist{29,tttttttt}));
    end 

    
    % Temperature dependence
    if ismissing(todolist{30, tttttttt}) || ismissing(todolist{31, tttttttt})
        assignin('base','TempDenp_tf','no');
    else 
        assignin('base','TempDenp_tf','yes');
    end   
    if ismissing(todolist{30, tttttttt})
        assignin('base','Temps_4_Denp',num2str([]));
    else
        assignin('base','Temps_4_Denp',num2str(todolist{30,tttttttt}));
    end 
    if ismissing(todolist{31, tttttttt})
        assignin('base','TempSensi_4_Denp',num2str([]));
    else
        assignin('base','TempSensi_4_Denp',num2str(todolist{31,tttttttt}));
    end 


    % tau and/or R relimit
    if ismissing(todolist{32, tttttttt}) && ismissing(todolist{33, tttttttt}) && ismissing(todolist{34, tttttttt}) && ismissing(todolist{35, tttttttt})
        assignin('base','Ri_taui_relimit_tf','no');
    else 
        assignin('base','Ri_taui_relimit_tf','yes');
    end   
    if ismissing(todolist{32, tttttttt})
        assignin('base','tau1_limit',num2str([]));
    else
        assignin('base','tau1_limit',num2str(todolist{32,tttttttt}));
    end 
    if ismissing(todolist{33, tttttttt})
        assignin('base','tau2_limit',num2str([]));
    else
        assignin('base','tau2_limit',num2str(todolist{33,tttttttt}));
    end 
    if ismissing(todolist{34, tttttttt})
        assignin('base','tau3_limit',num2str([]));
    else
        assignin('base','tau3_limit',num2str(todolist{34,tttttttt}));
    end 
    if ismissing(todolist{35, tttttttt})
        assignin('base','Ri_limit',num2str([]));
    else
        assignin('base','Ri_limit',num2str(todolist{35,tttttttt}));
    end 

    disp('###################################################################################################################')
    disp(['Datasheet ',num2str(tttttttt-2),' of ',num2str(column_dataend-2)])
    run ECM_Easy.m

    fig1h = findall(0,'type','figure','Tag','figure1');
    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, fig1h);
    delete(other_figures);  % close all figures except GUI
end

%% end of run
disp('###################################################################################################################')
disp('Congratulations!!! Batch parameterisation completed !!! Please go to the folder \Results\ to check all the results.')
disp('(For debugging (or donation) purpose, please contact tao.zhu@imperial.ac.uk :D)')

