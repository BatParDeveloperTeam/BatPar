cd Results
cd(fileNAME)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating LUTs with temperature dependence @ different Crates %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read-write LUT and generate plots
for iiii = 1 : length(C_RATE_DADDY)

    % unique SOC appeared in every temperature
    SOC_DADDY = [];
    for jjjj = 1: length(TEMP_DADDY)
        SOC_DADDY = [SOC_DADDY;LUT_DADDY{iiii,jjjj}(:,1)]; 
    end
    SOC_DADDY = sort(unique(SOC_DADDY));
    
    % Create LUTs
    LUT__R0{iiii} = ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
    LUT__R0{iiii}(2:end,1) = SOC_DADDY;
    LUT__R0{iiii}(1,2:end) = TEMP_DADDY;
    
    LUT__R1{iiii}= ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
    LUT__R1{iiii}(2:end,1) = SOC_DADDY;
    LUT__R1{iiii}(1,2:end) = TEMP_DADDY;
    
    LUT__C1{iiii}= ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
    LUT__C1{iiii}(2:end,1) = SOC_DADDY;
    LUT__C1{iiii}(1,2:end) = TEMP_DADDY;
    
    LUT__Tau1{iiii}= ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
    LUT__Tau1{iiii}(2:end,1) = SOC_DADDY;
    LUT__Tau1{iiii}(1,2:end) = TEMP_DADDY;
    
    if NRC == 2 || NRC == 3
        LUT__R2{iiii}= ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__R2{iiii}(2:end,1) = SOC_DADDY;
        LUT__R2{iiii}(1,2:end) = TEMP_DADDY;
        
        LUT__C2{iiii}= ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__C2{iiii}(2:end,1) = SOC_DADDY;
        LUT__C2{iiii}(1,2:end) = TEMP_DADDY;
        
        LUT__Tau2{iiii} = ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__Tau2{iiii}(2:end,1) = SOC_DADDY;
        LUT__Tau2{iiii}(1,2:end) = TEMP_DADDY;
    end
    
    if NRC == 3
        LUT__R3{iiii} = ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__R3{iiii}(2:end,1) = SOC_DADDY;
        LUT__R3{iiii}(1,2:end) = TEMP_DADDY;
        
        LUT__C3{iiii} = ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__C3{iiii}(2:end,1) = SOC_DADDY;
        LUT__C3{iiii}(1,2:end) = TEMP_DADDY;
        
        LUT__Tau3{iiii} = ones(length(SOC_DADDY)+1,length(TEMP_DADDY)+1)*nan;
        LUT__Tau3{iiii}(2:end,1) = SOC_DADDY;
        LUT__Tau3{iiii}(1,2:end) = TEMP_DADDY;
    end
    
    % Write LUTs
    for rrrr = 1 : length(TEMP_DADDY)   % number of pre-LUTs (temperature specific)
        for ssss = 2 : length(LUT__R0{iiii}(:,1))   % number of SOCs in post-LUTs
        index_soc = find (LUT_DADDY{iiii,rrrr}(:,1) == LUT__R0{iiii}(ssss,1));
            if ~isempty(index_soc)
                LUT__R0{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,3);
                LUT__R1{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,4);
                LUT__C1{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,7);
                LUT__Tau1{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,10);
                if NRC == 2 || NRC == 3
                    LUT__R2{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,5);
                    LUT__C2{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,8);
                    LUT__Tau2{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,11);
                end
                if NRC == 3
                    LUT__R3{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,6);
                    LUT__C3{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,9);
                    LUT__Tau3{iiii} (ssss , rrrr+1) = LUT_DADDY{iiii,rrrr}(index_soc,12);
                end 
            end
        end
    end 
    
    % export LUTs
    writematrix(LUT__R0{iiii},'LUT_R0_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__R1{iiii},'LUT_R1_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__C1{iiii},'LUT_C1_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__Tau1{iiii},'LUT_Tau1_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    
    if NRC == 2 || NRC == 3
    writematrix(LUT__R2{iiii},'LUT_R2_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__C2{iiii},'LUT_C2_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__Tau2{iiii},'LUT_Tau2_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    end
    
    if NRC == 3
    writematrix(LUT__R3{iiii},'LUT_R3_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__C3{iiii},'LUT_C3_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    writematrix(LUT__Tau3{iiii},'LUT_Tau3_with_Temperatures @ All C-rates.xlsx','Sheet',[num2str(C_RATE_DADDY(iiii)),'C']);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating plots with temperature dependence @ different Crates %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotrownumber = ceil(length(C_RATE_DADDY)/3);


figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(C_RATE_DADDY) 
    SOC_DADDY = LUT__R0{ttt}(2:end, 1)';

    [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
    R0_mesh = LUT__R0{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(TEMP_mesh, SOC_mesh , R0_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(TEMP_mesh, SOC_mesh , R0_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('R0 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['R0 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
end
savefig('Parameterised R0 with Temperatures @ different C-rates.fig');


figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(C_RATE_DADDY) 
    SOC_DADDY = LUT__R1{ttt}(2:end, 1)';

    [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
    R1_mesh = LUT__R1{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(TEMP_mesh, SOC_mesh , R1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(TEMP_mesh, SOC_mesh , R1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('R1 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['R1 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
end
savefig('Parameterised R1 with Temperatures @ different C-rates.fig');


figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(C_RATE_DADDY) 
    SOC_DADDY = LUT__C1{ttt}(2:end, 1)';

    [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
    C1_mesh = LUT__C1{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(TEMP_mesh, SOC_mesh , C1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(TEMP_mesh, SOC_mesh , C1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('C1 (Farad)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['C1 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
end
savefig('Parameterised C1 with Temperatures @ different C-rates.fig');


figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(C_RATE_DADDY) 
    SOC_DADDY = LUT__Tau1{ttt}(2:end, 1)';

    [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
    Tau1_mesh = LUT__Tau1{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(TEMP_mesh, SOC_mesh , Tau1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(TEMP_mesh, SOC_mesh , Tau1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('Tau1'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['Tau1 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
end
savefig('Parameterised Tau1 with Temperatures @ different C-rates.fig');

    
if NRC == 2 || NRC == 3
    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__R2{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        R2_mesh = LUT__R2{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , R2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , R2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('R2 (ohm)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['R2 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised R2 with Temperatures @ different C-rates.fig');
    

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__C2{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        C2_mesh = LUT__C2{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , C2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , C2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('C2 (Farad)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['C2 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised C2 with Temperatures @ different C-rates.fig');

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__Tau2{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        Tau2_mesh = LUT__Tau2{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , Tau2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , Tau2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('Tau2'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['Tau2 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised Tau2 with Temperatures @ different C-rates.fig');
end


if NRC == 3
    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__R3{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        R3_mesh = LUT__R3{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , R3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , R3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('R3 (ohm)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['R3 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised R3 with Temperatures @ different C-rates.fig');
    

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__C3{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        C3_mesh = LUT__C3{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , C3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , C3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('C3 (Farad)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['C3 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised C3 with Temperatures @ different C-rates.fig');

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(C_RATE_DADDY) 
        SOC_DADDY = LUT__Tau3{ttt}(2:end, 1)';
    
        [TEMP_mesh , SOC_mesh] = meshgrid(TEMP_DADDY, SOC_DADDY) ; % mesh it 
        Tau3_mesh = LUT__Tau3{ttt}(2:length(SOC_DADDY) +1 , 2:length(TEMP_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(TEMP_mesh, SOC_mesh , Tau3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(TEMP_mesh, SOC_mesh , Tau3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Temperature (degC)'); ylabel('SOC'); zlabel('Tau3'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['Tau3 with Temperatures @ ',num2str(C_RATE_DADDY(ttt)),'C']); 
    end
    savefig('Parameterised Tau3 with Temperatures @ different C-rates.fig');    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating LUTs with temperature dependence @ different Crates %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read-write LUT and generate plots
for iiii = 1 : length(TEMP_DADDY)

    % unique SOC appeared in every temperature
    SOC_DADDY = [];
    for jjjj = 1: length(C_RATE_DADDY)
        SOC_DADDY = [SOC_DADDY;LUT_DADDY{jjjj,iiii}(:,1)]; 
    end
    SOC_DADDY = sort(unique(SOC_DADDY));
    
    % Create LUTs
    LUT__R0{iiii} = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__R0{iiii}(2:end,1) = SOC_DADDY;
    LUT__R0{iiii}(1,2:end) = C_RATE_DADDY;
    
    LUT__R1{iiii}= ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__R1{iiii}(2:end,1) = SOC_DADDY;
    LUT__R1{iiii}(1,2:end) = C_RATE_DADDY;
    
    LUT__C1{iiii}= ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__C1{iiii}(2:end,1) = SOC_DADDY;
    LUT__C1{iiii}(1,2:end) = C_RATE_DADDY;
    
    LUT__Tau1{iiii}= ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__Tau1{iiii}(2:end,1) = SOC_DADDY;
    LUT__Tau1{iiii}(1,2:end) = C_RATE_DADDY;
    
    if NRC == 2 || NRC == 3
        LUT__R2{iiii}= ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__R2{iiii}(2:end,1) = SOC_DADDY;
        LUT__R2{iiii}(1,2:end) = C_RATE_DADDY;
        
        LUT__C2{iiii}= ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__C2{iiii}(2:end,1) = SOC_DADDY;
        LUT__C2{iiii}(1,2:end) = C_RATE_DADDY;
        
        LUT__Tau2{iiii} = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__Tau2{iiii}(2:end,1) = SOC_DADDY;
        LUT__Tau2{iiii}(1,2:end) = C_RATE_DADDY;
    end
    
    if NRC == 3
        LUT__R3{iiii} = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__R3{iiii}(2:end,1) = SOC_DADDY;
        LUT__R3{iiii}(1,2:end) = C_RATE_DADDY;
        
        LUT__C3{iiii} = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__C3{iiii}(2:end,1) = SOC_DADDY;
        LUT__C3{iiii}(1,2:end) = C_RATE_DADDY;
        
        LUT__Tau3{iiii} = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
        LUT__Tau3{iiii}(2:end,1) = SOC_DADDY;
        LUT__Tau3{iiii}(1,2:end) = C_RATE_DADDY;
    end
    
    % Write LUTs
    for rrrr = 1 : length(C_RATE_DADDY)   % number of pre-LUTs (current specific)
        for ssss = 2 : length(LUT__R0{iiii}(:,1))   % number of SOCs in post-LUTs
        index_soc = find (LUT_DADDY{rrrr,iiii}(:,1) == LUT__R0{iiii}(ssss,1));
            if ~isempty(index_soc)
                LUT__R0{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,3);
                LUT__R1{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,4);
                LUT__C1{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,7);
                LUT__Tau1{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,10);
                if NRC == 2 || NRC == 3
                    LUT__R2{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,5);
                    LUT__C2{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,8);
                    LUT__Tau2{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,11);
                end
                if NRC == 3
                    LUT__R3{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,6);
                    LUT__C3{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,9);
                    LUT__Tau3{iiii} (ssss , rrrr+1) = LUT_DADDY{rrrr,iiii}(index_soc,12);
                end 
            end
        end
    end 
    
    % Export LUTs
    writematrix(LUT__R0{iiii},'LUT_R0_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__R1{iiii},'LUT_R1_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__C1{iiii},'LUT_C1_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__Tau1{iiii},'LUT_Tau1_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    
    if NRC == 2 || NRC == 3
    writematrix(LUT__R2{iiii},'LUT_R2_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__C2{iiii},'LUT_C2_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__Tau2{iiii},'LUT_Tau2_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    end
    
    if NRC == 3
    writematrix(LUT__R3{iiii},'LUT_R3_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__C3{iiii},'LUT_C3_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    writematrix(LUT__Tau3{iiii},'LUT_Tau3_with_C-rates @ All Temperatures.xlsx','Sheet',[num2str(TEMP_DADDY(iiii)),'degC']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating plots with Crate dependence @ different temperatures %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotrownumber = ceil(length(TEMP_DADDY)/3);

figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(TEMP_DADDY) 
    SOC_DADDY = LUT__R0{ttt}(2:end, 1)';

    [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
    R0_mesh = LUT__R0{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(Crate_mesh, SOC_mesh , R0_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(Crate_mesh, SOC_mesh , R0_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current rate (C)'); ylabel('SOC'); zlabel('R0 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['R0 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
end
savefig('Parameterised R0 with C-rates @ different Temperatures.fig');


figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(TEMP_DADDY) 
    SOC_DADDY = LUT__R1{ttt}(2:end, 1)';

    [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
    R1_mesh = LUT__R1{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(Crate_mesh, SOC_mesh , R1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(Crate_mesh, SOC_mesh , R1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current rate (C)'); ylabel('SOC'); zlabel('R1 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['R1 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
end
savefig('Parameterised R1 with C-rates @ different Temperatures.fig');

figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(TEMP_DADDY) 
    SOC_DADDY = LUT__C1{ttt}(2:end, 1)';

    [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
    C1_mesh = LUT__C1{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(Crate_mesh, SOC_mesh , C1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(Crate_mesh, SOC_mesh , C1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current rate (C)'); ylabel('SOC'); zlabel('C1 (Farad)'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['C1 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
end
savefig('Parameterised C1 with C-rates @ different Temperatures.fig');

figure('Position', [200, 0, 1700, 450*plotrownumber]);
for ttt = 1: length(TEMP_DADDY) 
    SOC_DADDY = LUT__Tau1{ttt}(2:end, 1)';

    [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
    Tau1_mesh = LUT__Tau1{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    subplot (plotrownumber, 3 , ttt);
    grid on
    hold on

    hahaha = surf(Crate_mesh, SOC_mesh , Tau1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    scatter3(Crate_mesh, SOC_mesh , Tau1_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current rate (C)'); ylabel('SOC'); zlabel('Tau1'); 
    axis([-inf inf 0 1 -inf inf]);
    title (['Tau1 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
end
savefig('Parameterised Tau1 with C-rates @ different Temperatures.fig');

    
if NRC == 2 || NRC == 3
    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__R2{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        R2_mesh = LUT__R2{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , R2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , R2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('R2 (ohm)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['R2 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised R2 with C-rates @ different Temperatures.fig');
    

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__C2{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        C2_mesh = LUT__C2{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , C2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , C2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('C2 (Farad)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['C2 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised C2 with C-rates @ different Temperatures.fig');

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__Tau2{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        Tau2_mesh = LUT__Tau2{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , Tau2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , Tau2_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('Tau2'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['Tau2 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised Tau2 with C-rates @ different Temperatures.fig');
end

if NRC == 3
    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__R3{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        R3_mesh = LUT__R3{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , R3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , R3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('R3 (ohm)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['R3 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised R3 with C-rates @ different Temperatures.fig');
    

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__C3{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        C3_mesh = LUT__C3{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , C3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , C3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('C3 (Farad)'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['C3 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised C3 with C-rates @ different Temperatures.fig');

    figure('Position', [200, 0, 1700, 450*plotrownumber]);
    for ttt = 1: length(TEMP_DADDY) 
        SOC_DADDY = LUT__Tau3{ttt}(2:end, 1)';
    
        [Crate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY) ; % mesh it 
        Tau3_mesh = LUT__Tau3{ttt}(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
        
        subplot (plotrownumber, 3 , ttt);
        grid on
        hold on
    
        hahaha = surf(Crate_mesh, SOC_mesh , Tau3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
        hahaha.EdgeColor = 'none';
        scatter3(Crate_mesh, SOC_mesh , Tau3_mesh, 150, 'k.');
        view([-45 30]);
        xlabel('Current rate (C)'); ylabel('SOC'); zlabel('Tau3'); 
        axis([-inf inf 0 1 -inf inf]);
        title (['Tau3 with C-rates @ ',num2str(TEMP_DADDY(ttt)),'degC']); 
    end
    savefig('Parameterised Tau3 with C-rates @ different Temperatures.fig'); 
end