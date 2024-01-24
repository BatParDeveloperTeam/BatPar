function [dataIII, SOC_steps]= data_currentDependent (data, Currents_depended, Currents_sensitivity, S_j , NRC)

% data reorganisation
for i = 1:length(Currents_depended)  % every current
    for j = 1: length(data.I)        % every SOC
        count = 0;
        for k = 1: length(data.I{1,j}) % every segment
            
            indexIII = currentIndexing (data.I{1,j}{k}, Currents_depended(i), Currents_sensitivity);

            if any(indexIII)
                [~, minlocs] = findpeaks(-indexIII);
                if ~isempty(minlocs)
                    
                    count = count + 1;  
                    indexIII_temp = indexIII;
                    indexIII_temp(minlocs(1) : end) = false;
                    dataIII.I {i,j}{count,1} = data.I{1,j}{k}(indexIII_temp);
                    dataIII.t {i,j}{count,1} = data.t{1,j}{k}(indexIII_temp);
                    dataIII.T {i,j}{count,1} = data.T{1,j}{k}(indexIII_temp);
                    dataIII.V {i,j}{count,1} = data.V{1,j}{k}(indexIII_temp);
                    dataIII.S {i,j}{count,1} = data.S{1,j}{k}(indexIII_temp);
                    
                    count = count + 1;  
                    indexIII_temp = indexIII;
                    indexIII_temp(1 : minlocs(end)) = false;
                    dataIII.I {i,j}{count,1} = data.I{1,j}{k}(indexIII_temp);
                    dataIII.t {i,j}{count,1} = data.t{1,j}{k}(indexIII_temp);
                    dataIII.T {i,j}{count,1} = data.T{1,j}{k}(indexIII_temp);
                    dataIII.V {i,j}{count,1} = data.V{1,j}{k}(indexIII_temp);
                    dataIII.S {i,j}{count,1} = data.S{1,j}{k}(indexIII_temp);
                    
                    if  length(minlocs) > 1
                        for hhh = 2: length(minlocs)
                            count = count + 1;
                            indexIII_temp = indexIII;
                            indexIII_temp(1 : minlocs(hhh - 1)) = false;
                            indexIII_temp(minlocs(hhh) : end) = false;
                            dataIII.I {i,j}{count,1} = data.I{1,j}{k}(indexIII_temp);
                            dataIII.t {i,j}{count,1} = data.t{1,j}{k}(indexIII_temp);
                            dataIII.T {i,j}{count,1} = data.T{1,j}{k}(indexIII_temp);
                            dataIII.V {i,j}{count,1} = data.V{1,j}{k}(indexIII_temp);
                            dataIII.S {i,j}{count,1} = data.S{1,j}{k}(indexIII_temp);
                        end 
                    end 

                else
                    count = count + 1;   
                    dataIII.I {i,j}{count,1} = data.I{1,j}{k}(indexIII);
                    dataIII.t {i,j}{count,1} = data.t{1,j}{k}(indexIII);
                    dataIII.T {i,j}{count,1} = data.T{1,j}{k}(indexIII);
                    dataIII.V {i,j}{count,1} = data.V{1,j}{k}(indexIII);
                    dataIII.S {i,j}{count,1} = data.S{1,j}{k}(indexIII);
                end 
            end 
        end
    end
end


% data check --> filter out empty data cells

indexDataI = ~cellfun(@isempty,dataIII.t);

for zhu = 1:length(Currents_depended)
                
    dataIII_t = dataIII.t(zhu,:);
    dataIII_temp.t {zhu,1} = {dataIII_t(indexDataI(zhu,:))};         
    dataIII_I = dataIII.I(zhu,:);
    dataIII_temp.I {zhu,1} = {dataIII_I(indexDataI(zhu,:))};
    dataIII_T = dataIII.T(zhu,:);
    dataIII_temp.T {zhu,1} = {dataIII_T(indexDataI(zhu,:))};
    dataIII_V = dataIII.V(zhu,:);
    dataIII_temp.V {zhu,1} = {dataIII_V(indexDataI(zhu,:))};
    dataIII_S = dataIII.S(zhu,:);
    dataIII_temp.S {zhu,1} = {dataIII_S(indexDataI(zhu,:))};
    dataIII_temp.S_j_I {zhu,1} = {S_j(indexDataI(zhu,:))};
        
end

% clean down
for tao = 1:length(Currents_depended)
    for i = 1: length(dataIII_temp.t{tao,1}{1,1})
        ind = cellfun('length', dataIII_temp.t{tao,1}{1,1}{1,i}) >= (NRC * 3 + 2) * 2;
        dataIII_temp_temp.t{tao,1}{1,1}{1,i} = dataIII_temp.t{tao,1}{1,1}{1,i}(ind);
        dataIII_temp_temp.I{tao,1}{1,1}{1,i} = dataIII_temp.I{tao,1}{1,1}{1,i}(ind);
        dataIII_temp_temp.T{tao,1}{1,1}{1,i} = dataIII_temp.T{tao,1}{1,1}{1,i}(ind);
        dataIII_temp_temp.V{tao,1}{1,1}{1,i} = dataIII_temp.V{tao,1}{1,1}{1,i}(ind);
        dataIII_temp_temp.S{tao,1}{1,1}{1,i} = dataIII_temp.S{tao,1}{1,1}{1,i}(ind);
    end
end 
dataIII_temp_temp.S_j_I = dataIII_temp.S_j_I;

dataIII = dataIII_temp_temp;

% steps for each current
SOC_steps = [];
for jj = 1: length(Currents_depended)
   SOC_steps (jj,1) = length(dataIII.S_j_I{jj,1}{1,1});
end

end