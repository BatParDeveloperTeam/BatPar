function data = cleandown(data,NRC)

for iii = 1:numel(data.t)
    
    ind = cellfun('length',data.t{iii})>= (NRC * 3 + 2) * 2; % cellfun - apply func 'length' to !!!!! each matrix of each cell !!!!! in data.t. Return a cell of same dimension
    % Why (NRC * 3 + 2) * 2? Because:
    % Firstly, there are (NRC * 3 + 2) unknowns in the ECM, as each RC pair has three unknowns - Ri, tau_i and V0_i, and the ECM has two more unknowns - OCV and R0
    % Secondly, * 2 is a safety coefficient to ensure a unique solution. i.e., the solver (fmincon) is fed with one more fold data than necessary
  
    data.t{iii} = data.t{iii}(ind);   % if ind [ijk]= 1, keep data.t {iii}[ijk]; if ind [ijk]=0, remove data.t {iii}[ijk] 
    data.I{iii} = data.I{iii}(ind);
    data.V{iii} = data.V{iii}(ind);    
    data.T{iii} = data.T{iii}(ind);
    data.S{iii} = data.S{iii}(ind);
    data.label{iii} = data.label{iii}(ind);
    
end

end