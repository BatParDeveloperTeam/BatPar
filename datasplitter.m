function [data,S_j] = datasplitter(t,I,V,T,S,S_j,bounds,flags)

OPT_factor = 1;
% 
% if (strcmp(flags.SOCOCVtable.tf,'yes'))  % string compare
%      OPT_factor = 2;   
% end

% Lower and upper bounds of each SOC window;
for ii = 1:numel(S_j)
    
        if strcmp(flags.MaxSOCis1.tf,'yes')  % 100% SOC reached
     
           if (ii-1==0)
               S_pesudobound = S_j(ii) - (S_j(ii+2) - S_j(ii+1))/OPT_factor;
               if bounds.lowerbound > S_pesudobound
		        lowerbound = bounds.lowerbound;
		       else
		        lowerbound = S_pesudobound;
		       end
           else 
               lowerbound = S_j(ii) - (S_j(ii) - S_j(ii-1))/OPT_factor;
           end
           
           if (ii+1>numel(S_j))
               upperbound = S_j(ii);                           % is 1 (if the cell was ever fully charged)
               % lowerbound = S_j(ii) - (S_j(ii) - S_j(ii-1));
           else
               upperbound = S_j(ii) + (S_j(ii+1)-S_j(ii))/OPT_factor;   
           end
        
        else  % 0% SOC reached
           
           if (ii-1==0)          
		       lowerbound = S_j(ii);                          % is 0 (if the cell was ever fully discharged)
           else 
               lowerbound = S_j(ii) - (S_j(ii) - S_j(ii-1))/OPT_factor;
           end
           
           if (ii+1>numel(S_j))
               S_pesudobound = S_j(ii) + (S_j(ii-1) - S_j(ii-2))/OPT_factor;
               if bounds.upperbound < S_pesudobound
		        upperbound = bounds.upperbound;
		       else
		        upperbound = S_pesudobound;
		       end                          
           else
               upperbound = S_j(ii) + (S_j(ii+1)-S_j(ii))/OPT_factor;   
           end
        end 
       % Find indices to find data in the given window
  
       count = 0;
       datatemp = [];  % Fixed on 16 June 2022 --> wipe off previous calculation !

       for ijk = 1:numel(t) % for each charge segment

           ind1 = (lowerbound <= S{ijk}) + (S{ijk} <= upperbound) == 2;     % every element in S{} will be compared with lowerbound and upperbound
           % (lowerbound <= S{ijk}) + (S{ijk} <= upperbound) -- a vector ranges in (0,1,2). Element value in this vector is 2 when both conditions satisfied.
           % Whihch means the element lies in between lower and upper bound
           % ind1 logical (zero-one) vector
           if max(ind1)>0 % if there is at least one element lying in this SOC segment
  
                ind2 = cumsum(ones(size(ind1)));  % 1:1:lengh(ind)
                ind2 = ind2(ind1);   % keep the line number of S{} where elements lie in bounds 

                datatemp.t{1+count,1} = t{ijk}(ind1);
                datatemp.I{1+count,1} = I{ijk}(ind1);
                datatemp.V{1+count,1} = V{ijk}(ind1); % mapping ind2 to t, I, V, T, S
                datatemp.T{1+count,1} = T{ijk}(ind1);
                datatemp.S{1+count,1} = S{ijk}(ind1);
                datatemp.label{1+count,1} = [ijk,ind2(1),ind2(end)]; % [segment number, start index, end index]
  
                data.t{ii}        = datatemp.t;
                data.I{ii}        = datatemp.I;
                data.V{ii}        = datatemp.V;
                data.T{ii}        = datatemp.T;
                data.S{ii}        = datatemp.S;
                data.label{ii}    = datatemp.label;

                count = count + 1;
            else
		data.t{ii}        = [];
                data.I{ii}        = [];
                data.V{ii}        = [];
                data.T{ii}        = [];
                data.S{ii}        = [];
                data.label{ii}    = [];
                
                count = count + 1;  	
            end
 
        end


end

% data check --> filter out empty data clusters
count = 0;
datainsufficient = zeros(1,numel(S_j));
for zhu = 1:numel(S_j)
    if isempty (data.t{zhu})
        datainsufficient (zhu) = 1;
        count = count + 1;
    end
end
if count ~= 0
    datainsufficient = logical(datainsufficient);
%     disp('*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******')
%     disp(['Parameterisation cannot happen at SOC = ',num2str(S_j(datainsufficient)),' because of insufficient experimental data. The codes have to skip the mentioned SOC. If you keep seeing this message, please make the SOC window size larger.'])
%     disp('*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******CAUTION*******')
    
    S_j           = S_j(~datainsufficient);
    data.t        = data.t(~datainsufficient);
    data.I        = data.I(~datainsufficient);
    data.V        = data.V(~datainsufficient);
    data.T        = data.T(~datainsufficient);
    data.S        = data.S(~datainsufficient);
    data.label    = data.label(~datainsufficient);
end

end
  
