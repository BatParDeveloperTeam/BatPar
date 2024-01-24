function [dataTTT, SOC_steps]= data_temperatureDependent (data, Temps_depended, Temps_sensitivity, S_j, NRC)

% data reorganisation
for i = 1:length(Temps_depended)     % every Temp
    for j = 1: length(data.T)        % every SOC
        count = 0;
        for k = 1: length(data.T{1,j}) % every segment
            
            indexTTT =  ((Temps_depended(i) - Temps_sensitivity) <= data.T{1,j}{k}) + (data.T{1,j}{k} <= (Temps_depended(i) + Temps_sensitivity)) == 2;

            if any(indexTTT)
                    count = count + 1;   
                    dataTTT.I {i,j}{count,1} = data.I{1,j}{k}(indexTTT);
                    dataTTT.t {i,j}{count,1} = data.t{1,j}{k}(indexTTT);
                    dataTTT.T {i,j}{count,1} = data.T{1,j}{k}(indexTTT);
                    dataTTT.V {i,j}{count,1} = data.V{1,j}{k}(indexTTT);
                    dataTTT.S {i,j}{count,1} = data.S{1,j}{k}(indexTTT);
            end 
        end
    end
end

% data check --> filter out empty data cells
indexDataT = ~cellfun(@isempty,dataTTT.t);

for zhu = 1:length(Temps_depended)
                
    dataTTT_t = dataTTT.t(zhu,:);
    dataTTT_temp.t {zhu,1} = {dataTTT_t(indexDataT(zhu,:))};         
    dataTTT_I = dataTTT.I(zhu,:);
    dataTTT_temp.I {zhu,1} = {dataTTT_I(indexDataT(zhu,:))};
    dataTTT_T = dataTTT.T(zhu,:);
    dataTTT_temp.T {zhu,1} = {dataTTT_T(indexDataT(zhu,:))};
    dataTTT_V = dataTTT.V(zhu,:);
    dataTTT_temp.V {zhu,1} = {dataTTT_V(indexDataT(zhu,:))};
    dataTTT_S = dataTTT.S(zhu,:);
    dataTTT_temp.S {zhu,1} = {dataTTT_S(indexDataT(zhu,:))};
    dataTTT_temp.S_j_I {zhu,1} = {S_j(indexDataT(zhu,:))};
        
end

% clean down
for tao = 1:length(Temps_depended)
    for i = 1: length(dataTTT_temp.t{tao,1}{1,1})
        ind = cellfun('length', dataTTT_temp.t{tao,1}{1,1}{1,i}) >= (NRC * 3 + 2) * 2;
        dataTTT_temp_temp.t{tao,1}{1,1}{1,i} = dataTTT_temp.t{tao,1}{1,1}{1,i}(ind);
        dataTTT_temp_temp.I{tao,1}{1,1}{1,i} = dataTTT_temp.I{tao,1}{1,1}{1,i}(ind);
        dataTTT_temp_temp.T{tao,1}{1,1}{1,i} = dataTTT_temp.T{tao,1}{1,1}{1,i}(ind);
        dataTTT_temp_temp.V{tao,1}{1,1}{1,i} = dataTTT_temp.V{tao,1}{1,1}{1,i}(ind);
        dataTTT_temp_temp.S{tao,1}{1,1}{1,i} = dataTTT_temp.S{tao,1}{1,1}{1,i}(ind);
    end
end 
dataTTT_temp_temp.S_j_I = dataTTT_temp.S_j_I;

dataTTT = dataTTT_temp_temp;

% steps for each current
SOC_steps = [];
for jj = 1: length(Temps_depended)
   SOC_steps (jj,1) = length(dataTTT.S_j_I{jj,1}{1,1});
end

end
