clc
cd(evalin('base','exp_path'));

Filescsv = dir(strcat('*.csv'));
Filesxlsx = dir(strcat('*.xlsx'));

Files = [Filesxlsx;Filescsv];
disp(['Batch run started, using inputs defined in the GUI. ',num2str(length(Files)),' datasheets found. Good luck!'])
assignin('base','Files',Files);

cd ..


for tttttttt = 1:length(Files)
    clearvars -except tttttttt Files exp_path record_column time_column current_column voltage_column temp_column row_start row_end Nom_capacity voltage_up_limit voltage_low_limit OCVsheet_name RCpair_number SOC_window SOC_param_lowlimit SOC_param_uplimit tau1_limit tau2_limit tau3_limit Ri_limit Currents_4_Denp CurrentSensi_4_Denp Temps_4_Denp TempSensi_4_Denp Current_belwozero_tf MaxSOCis1 Mode_anode_tf Mode_charge_tf Mode_cathode_tf Pulse_head_tf fixedtau_realy_tf OCVtable_tf datasheet_cleanup_tf Ri_taui_relimit_tf CurrentDenp_tf TempDenp_tf 
    files_ = evalin('base','Files');
    disp('###################################################################################################################')
    disp(['Datasheet ',num2str(tttttttt),' of ',num2str(length(files_))])
    file_name = files_(tttttttt).name; 
    assignin('base','datasheet_name',file_name);
    run ECM_Easy.m
    
    fig1h = findall(0,'type','figure','Tag','figure1');
    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, fig1h);
    delete(other_figures);  % close all figures except GUI
end

disp('#######################################################################################################################')
disp('Congratulations!!! Batch parameterisation completed !!! Please go to the folder \Results\ to check all the results.')
disp('(For debugging (or donation) purpose, please contact tao.zhu@imperial.ac.uk :D)')
