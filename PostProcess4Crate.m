cd Results
cd(fileNAME)

% unique SOC appeared in every C-rate
SOC_DADDY = [];
for gggg = 1 : length(C_RATE_DADDY)
    SOC_DADDY = [SOC_DADDY;LUT_DADDY{gggg}(:,1)]; 
end
SOC_DADDY = sort(unique(SOC_DADDY));

% Create LUTs
LUT__R0 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
LUT__R0(2:end,1) = SOC_DADDY;
LUT__R0(1,2:end) = C_RATE_DADDY;

LUT__R1 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
LUT__R1(2:end,1) = SOC_DADDY;
LUT__R1(1,2:end) = C_RATE_DADDY;

LUT__C1 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
LUT__C1(2:end,1) = SOC_DADDY;
LUT__C1(1,2:end) = C_RATE_DADDY;

LUT__Tau1 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
LUT__Tau1(2:end,1) = SOC_DADDY;
LUT__Tau1(1,2:end) = C_RATE_DADDY;

if NRC == 2 || NRC == 3
    LUT__R2 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__R2(2:end,1) = SOC_DADDY;
    LUT__R2(1,2:end) = C_RATE_DADDY;
    
    LUT__C2 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__C2(2:end,1) = SOC_DADDY;
    LUT__C2(1,2:end) = C_RATE_DADDY;
    
    LUT__Tau2 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__Tau2(2:end,1) = SOC_DADDY;
    LUT__Tau2(1,2:end) = C_RATE_DADDY;
end

if NRC == 3
    LUT__R3 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__R3(2:end,1) = SOC_DADDY;
    LUT__R3(1,2:end) = C_RATE_DADDY;
    
    LUT__C3 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__C3(2:end,1) = SOC_DADDY;
    LUT__C3(1,2:end) = C_RATE_DADDY;
    
    LUT__Tau3 = ones(length(SOC_DADDY)+1,length(C_RATE_DADDY)+1)*nan;
    LUT__Tau3(2:end,1) = SOC_DADDY;
    LUT__Tau3(1,2:end) = C_RATE_DADDY;
end

% Write LUTs
for rrrr = 1 : length(LUT_DADDY)   % number of pre-LUTs
    for ssss = 2 : length(LUT__R0(:,1))   % number of SOCs in post-LUTs
    index_soc = find (LUT_DADDY{rrrr}(:,1) == LUT__R0(ssss,1));
        if ~isempty(index_soc)
            LUT__R0 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,3);
            LUT__R1 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,4);
            LUT__C1 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,7);
            LUT__Tau1 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,10);
            if NRC == 2 || NRC == 3
                LUT__R2 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,5);
                LUT__C2 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,8);
                LUT__Tau2 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,11);
            end
            if NRC == 3
                LUT__R3 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,6);
                LUT__C3 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,9);
                LUT__Tau3 (ssss , rrrr+1) = LUT_DADDY{rrrr}(index_soc,12);
            end 
        end
    end
end 

% export LUTs
writematrix(LUT__R0,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','R0');
writematrix(LUT__R1,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','R1');
writematrix(LUT__C1,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','C1');
writematrix(LUT__Tau1,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','Tau1');

if NRC == 2 || NRC == 3
writematrix(LUT__R2,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','R2');
writematrix(LUT__C2,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','C2');
writematrix(LUT__Tau2,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','Tau2');
end

if NRC == 3
writematrix(LUT__R3,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','R3');
writematrix(LUT__C3,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','C3');
writematrix(LUT__Tau3,'LookUpTable_All Params @ All C_rates.xlsx','Sheet','Tau3');
end

% plots
[C_rate_mesh , SOC_mesh] = meshgrid(C_RATE_DADDY, SOC_DADDY') ; %mesh it 
R0_mesh = LUT__R0(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
R1_mesh = LUT__R1(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
C1_mesh = LUT__C1(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
Tau1_mesh = LUT__Tau1(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);

figure; grid on
hahaha = surf(C_rate_mesh, SOC_mesh , R0_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
hahaha.EdgeColor = 'none';
hold on
scatter3(C_rate_mesh, SOC_mesh , R0_mesh, 150, 'k.');
view([-45 30]);
xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('R0 (ohm)'); 
axis([-inf inf 0 1 -inf inf]);
title ('Parameterised R0 with SOC & C-rate'); 
savefig('Parameterised R0 @ All C-rates.fig');

figure; grid on
hahaha = surf(C_rate_mesh, SOC_mesh , R1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
hahaha.EdgeColor = 'none';
hold on
scatter3(C_rate_mesh, SOC_mesh , R1_mesh, 150, 'k.');
view([-45 30]);
xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('R1 (ohm)'); 
axis([-inf inf 0 1 -inf inf]);
title ('Parameterised R1 with SOC & C-rate'); 
savefig('Parameterised R1 @ All C-rates.fig');

figure; grid on
hahaha = surf(C_rate_mesh, SOC_mesh , C1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
hahaha.EdgeColor = 'none';
hold on
scatter3(C_rate_mesh, SOC_mesh , C1_mesh, 150, 'k.');
view([-45 30]);
xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('C1 (Farad)'); 
axis([-inf inf 0 1 -inf inf]);
title ('Parameterised C1 with SOC & C-rate'); 
savefig('Parameterised C1 @ All C-rates.fig');

figure; grid on
hahaha = surf(C_rate_mesh, SOC_mesh , Tau1_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
hahaha.EdgeColor = 'none';
hold on
scatter3(C_rate_mesh, SOC_mesh , Tau1_mesh, 150, 'k.');
view([-45 30]);
xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('Tau1'); 
axis([-inf inf 0 1 -inf inf]);
title ('Parameterised Tau1 with SOC & C-rate'); 
savefig('Parameterised Tau1 @ All C-rates.fig');

if NRC == 2 || NRC == 3
    R2_mesh = LUT__R2(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    C2_mesh = LUT__C2(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    Tau2_mesh = LUT__Tau2(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , R2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    hold on
    scatter3(C_rate_mesh, SOC_mesh , R2_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('R2 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised R2 with SOC & C-rate'); 
    savefig('Parameterised R2 @ All C-rates.fig');
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , C2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    hold on
    scatter3(C_rate_mesh, SOC_mesh , C2_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('C2 (Farad)'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised C2 with SOC & C-rate'); 
    savefig('Parameterised C2 @ All C-rates.fig');
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , Tau2_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    hold on
    scatter3(C_rate_mesh, SOC_mesh , Tau2_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('Tau2'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised Tau2 with SOC & C-rate'); 
    savefig('Parameterised Tau2 @ All C-rates.fig');
end

if NRC == 3
    R3_mesh = LUT__R3(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    C3_mesh = LUT__C3(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    Tau3_mesh = LUT__Tau3(2:length(SOC_DADDY) +1 , 2:length(C_RATE_DADDY)+1);
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , R3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    hold on
    scatter3(C_rate_mesh, SOC_mesh , R3_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('R3 (ohm)'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised R3 with SOC & C-rate'); 
    savefig('Parameterised R3 @ All C-rates.fig');
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , C3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    view([-45 30]);
    hold on
    scatter3(C_rate_mesh, SOC_mesh , C3_mesh, 150, 'k.');
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('C3 (Farad)'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised C3 with SOC & C-rate'); 
    savefig('Parameterised C3 @ All C-rates.fig');
    
    figure; grid on
    hahaha = surf(C_rate_mesh, SOC_mesh , Tau3_mesh,'FaceAlpha',0.95,'FaceColor','interp');        % fitting surface 
    hahaha.EdgeColor = 'none';
    hold on
    scatter3(C_rate_mesh, SOC_mesh , Tau3_mesh, 150, 'k.');
    view([-45 30]);
    xlabel('Current Rate (C)'); ylabel('SOC'); zlabel('Tau3'); 
    axis([-inf inf 0 1 -inf inf]);
    title ('Parameterised Tau3 with SOC & C-rate'); 
    savefig('Parameterised Tau3 @ All C-rates.fig');
end