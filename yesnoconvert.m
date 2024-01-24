function flags = yesnoconvert(flags)

    aux = fields(flags);  % get the attributes of flags; return a vector of attribute names
    
    for kk = 1:numel(aux)
        if (strcmp(flags.(aux{kk}).tf,'yes'))  % string compare
            flags.(aux{kk}).tf = 1;
        else
            flags.(aux{kk}).tf = 0;
        end
    end

end

