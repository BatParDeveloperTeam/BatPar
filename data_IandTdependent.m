function [dataIT, SOC_steps]= data_IandTdependent (data, Currents_depended, Currents_sensitivity, Temps_depended, Temps_sensitivity, S_j , NRC)

% data reorganisation - specific to I
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

% data reorganisation - specific to T
for ttt = 1: length(Temps_depended)            % every temperature 
    for i = 1:length(Currents_depended)      % every current
      for j = 1:length(S_j)                  % every SOC  
          count = 0;
          for k = 1:length(dataIII.T{i,j})   
            index_IT = ((Temps_depended(ttt) - Temps_sensitivity) <= dataIII.T{i,j}{k,1}) + (dataIII.T{i,j}{k,1} <= (Temps_depended(ttt) + Temps_sensitivity)) == 2;
            if any(index_IT) && (sum(index_IT) >= (NRC * 3 + 2) * 2)   % && clean down
                count = count + 1; 
                data_temp.I{i,ttt}{1,j}{count,1} = dataIII.I{i,j}{k,1}(index_IT);
                data_temp.t{i,ttt}{1,j}{count,1} = dataIII.t{i,j}{k,1}(index_IT);
                data_temp.T{i,ttt}{1,j}{count,1} = dataIII.T{i,j}{k,1}(index_IT);
                data_temp.V{i,ttt}{1,j}{count,1} = dataIII.V{i,j}{k,1}(index_IT);
                data_temp.S{i,ttt}{1,j}{count,1} = dataIII.S{i,j}{k,1}(index_IT);
            end
          end
          if count == 0
                data_temp.I{i,ttt}{1,j} = [];
                data_temp.t{i,ttt}{1,j} = [];
                data_temp.T{i,ttt}{1,j} = [];
                data_temp.V{i,ttt}{1,j} = [];
                data_temp.S{i,ttt}{1,j} = [];
          end
      end 
    end 
end 


% data check --> filter out empty data cells, plus number of steps for each IT
SOC_steps = [];
for ii = 1: length(Currents_depended)
    for jj = 1: length(Temps_depended)
        indexDataIT = ~cellfun(@isempty,data_temp.t{ii,jj});
        dataIT_temp.t {ii,jj} = data_temp.t{ii,jj}(indexDataIT);
        dataIT_temp.I {ii,jj} = data_temp.I{ii,jj}(indexDataIT);
        dataIT_temp.V {ii,jj} = data_temp.V{ii,jj}(indexDataIT);
        dataIT_temp.S {ii,jj} = data_temp.S{ii,jj}(indexDataIT);
        dataIT_temp.T {ii,jj} = data_temp.T{ii,jj}(indexDataIT);
        dataIT_temp.S_j_IT {ii,jj} = S_j(indexDataIT);
        SOC_steps (ii,jj) = length(dataIT_temp.S_j_IT{ii,jj});
    end
end

dataIT = dataIT_temp;

end