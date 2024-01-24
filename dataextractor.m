function [t, I, V, S, T, t_charge, I_charge, V_charge, S_charge, T_charge] =  dataextractor (t_, I_, V_, S_, T_)
    
%% Find local peaks as breakpoints
[~,maxlocs] = findpeaks(S_); % find values and locations of SOC local maximums (rising edge)
[minvals,minlocs] = findpeaks(-S_); 
minvals = -minvals;  % find values and locations of SOC local minimums (dropping edge)

if isempty(maxlocs) ==1 && isempty(minlocs) ==1  % single segment of charge or discharge 
    if mean(I_) > 0  % dischagre
       for pulse_head = 1:length(S_)-1  % find pulse_head - where the discharge pulse started
        diffS = S_(pulse_head+1)-S_(pulse_head);
        if diffS ~= 0
           break
        end 
       end

        t{1,1} = t_(pulse_head:end);
        I{1,1} = I_(pulse_head:end);
        V{1,1} = V_(pulse_head:end);
        T{1,1} = T_(pulse_head:end);
        S{1,1} = S_(pulse_head:end); 

        t_charge{1,1} = nan;
        I_charge{1,1} = nan;
        V_charge{1,1} = nan;
        T_charge{1,1} = nan;
        S_charge{1,1} = nan;
    else             % chagre
        for pulse_head = 1:length(S_)-1  % find pulse_head - where the charge pulse started
        diffS = S_(pulse_head+1)-S_(pulse_head);
            if diffS ~= 0            
            break
            end 
        end
        
        t{1,1} = nan;
        I{1,1} = nan;
        V{1,1} = nan;
        T{1,1} = nan;
        S{1,1} = nan; 

        t_charge{1,1} = t_(pulse_head:end);
        I_charge{1,1} = I_(pulse_head:end);
        V_charge{1,1} = V_(pulse_head:end);
        T_charge{1,1} = T_(pulse_head:end);
        S_charge{1,1} = S_(pulse_head:end);
    end
end

if isempty(maxlocs) ==0 && isempty(minlocs) ==1 % first chagre segment then discharge segment
    % shift the location of SOC local maximum from rising edge to dropping edge
    for nn = maxlocs(1):length(S_)-1
       diffS = S_(nn+1)-S_(maxlocs(1));
       if diffS ~= 0
           maxlocs(1) = nn;
           break
       end 
    end
    
    t{1,1} = t_(maxlocs(1):length(S_));
    I{1,1} = I_(maxlocs(1):length(S_));
    V{1,1} = V_(maxlocs(1):length(S_));
    T{1,1} = T_(maxlocs(1):length(S_));
    S{1,1} = S_(maxlocs(1):length(S_)); 
    

    for pulse_head = 1:maxlocs(1)  % find pulse_head - where the charge pulse started
        diffS = S_(pulse_head+1)-S_(pulse_head);
            if diffS ~= 0            
            break
            end 
    end
    t_charge{1,1} = t_(pulse_head:maxlocs(1));
    I_charge{1,1} = I_(pulse_head:maxlocs(1));
    V_charge{1,1} = V_(pulse_head:maxlocs(1));
    T_charge{1,1} = T_(pulse_head:maxlocs(1));
    S_charge{1,1} = S_(pulse_head:maxlocs(1));
end

if isempty(maxlocs) ==1 && isempty(minlocs) ==0  % first dischagre segment then charge segment
    % shift the location of SOC local minimum from dropping edge to rising edge
    for nn = minlocs(1):length(S_)-1
       diffS = S_(nn+1)-S_(minlocs(1));
       if diffS ~= 0
           minlocs(1) = nn;
           break
       end 
    end
    
    for pulse_head = 1:length(S_)-1  % find pulse_head - where the discharge pulse started
        diffS = S_(pulse_head+1)-S_(pulse_head);
            if diffS ~= 0            
            break
            end 
    end
    t{1,1} = t_(pulse_head:minlocs(1));
    I{1,1} = I_(pulse_head:minlocs(1));
    V{1,1} = V_(pulse_head:minlocs(1));
    T{1,1} = T_(pulse_head:minlocs(1));
    S{1,1} = S_(pulse_head:minlocs(1)); 
    
    t_charge{1,1} = t_(minlocs(1):length(S_));
    I_charge{1,1} = I_(minlocs(1):length(S_));
    V_charge{1,1} = V_(minlocs(1):length(S_));
    T_charge{1,1} = T_(minlocs(1):length(S_));
    S_charge{1,1} = S_(minlocs(1):length(S_));
end

if isempty(maxlocs) ==0 && isempty(minlocs) ==0 % multiple segments of both discharge and chagre     
for mm = 1:length(maxlocs) % shift the locations of SOC local maximums from rising edge to dropping edge
    for nn = maxlocs(mm):length(S_)-1
       diffS = S_(nn+1)-S_(maxlocs(mm));
       if diffS ~= 0
           maxlocs(mm) = nn;
           break
       end 
    end
end
       
for mm = 1:length(minlocs) % shift the locations of SOC local minimums from dropping edge to rising edge
    for nn = minlocs(mm):length(S_)-1
       diffS = S_(nn+1)-S_(minlocs(mm));
       if diffS ~= 0
           minlocs(mm) = nn;
           break
       end 
    end
end

%% identify wavefrom - number and order of peaks; split dicharge and charge segments

whocomesfirst=maxlocs(1)-minlocs(1); % local maximum or local minimum comes first?

if whocomesfirst < 0 % local maximum comes first
    if length(maxlocs)==length(minlocs) % same number
        for kk = 1:length(maxlocs)
            t{kk,1} = t_(maxlocs(kk):minlocs(kk));
            I{kk,1} = I_(maxlocs(kk):minlocs(kk));
            V{kk,1} = V_(maxlocs(kk):minlocs(kk));
            T{kk,1} = T_(maxlocs(kk):minlocs(kk));
            S{kk,1} = S_(maxlocs(kk):minlocs(kk)); % dischagre segments
            
            if kk == 1                             % chagre segments
                for pulse_head = 1:maxlocs(1)  % find pulse_head - where the charge pulse started
                    diffS = S_(pulse_head+1)-S_(pulse_head);
                    if diffS ~= 0            
                    break
                    end 
                end
                t_charge{kk,1} = t_(pulse_head:maxlocs(kk));
                I_charge{kk,1} = I_(pulse_head:maxlocs(kk));
                V_charge{kk,1} = V_(pulse_head:maxlocs(kk));
                T_charge{kk,1} = T_(pulse_head:maxlocs(kk));
                S_charge{kk,1} = S_(pulse_head:maxlocs(kk));
            else 
                if kk == length(maxlocs)
                    t_charge{kk+1,1} = t_(minlocs(kk):length(S_));
                    I_charge{kk+1,1} = I_(minlocs(kk):length(S_));
                    V_charge{kk+1,1} = V_(minlocs(kk):length(S_));
                    T_charge{kk+1,1} = T_(minlocs(kk):length(S_));
                    S_charge{kk+1,1} = S_(minlocs(kk):length(S_));
                else 
                    t_charge{kk,1} = t_(minlocs(kk-1):maxlocs(kk));
                    I_charge{kk,1} = I_(minlocs(kk-1):maxlocs(kk));
                    V_charge{kk,1} = V_(minlocs(kk-1):maxlocs(kk));
                    T_charge{kk,1} = T_(minlocs(kk-1):maxlocs(kk));
                    S_charge{kk,1} = S_(minlocs(kk-1):maxlocs(kk));
                end 
            end
        end
    else % number of local maximums (one) more than that of local minimums
        for kk = 1:length(maxlocs)-1
            t{kk,1} = t_(maxlocs(kk):minlocs(kk));
            I{kk,1} = I_(maxlocs(kk):minlocs(kk));
            V{kk,1} = V_(maxlocs(kk):minlocs(kk));
            T{kk,1} = T_(maxlocs(kk):minlocs(kk));
            S{kk,1} = S_(maxlocs(kk):minlocs(kk));
        end 
        t{length(maxlocs),1} = t_(maxlocs(length(maxlocs)):length(S_));
        I{length(maxlocs),1} = I_(maxlocs(length(maxlocs)):length(S_));
        V{length(maxlocs),1} = V_(maxlocs(length(maxlocs)):length(S_));
        T{length(maxlocs),1} = T_(maxlocs(length(maxlocs)):length(S_));
        S{length(maxlocs),1} = S_(maxlocs(length(maxlocs)):length(S_)); % dischagre segments
        
        for kk = 1:length(maxlocs)                                      % dischagre segments
            if kk == 1
                for pulse_head = 1:maxlocs(1)  % find pulse_head - where the charge pulse started
                    diffS = S_(pulse_head+1)-S_(pulse_head);
                    if diffS ~= 0            
                    break
                    end 
                end
                t_charge{kk,1} = t_(pulse_head:maxlocs(kk));
                I_charge{kk,1} = I_(pulse_head:maxlocs(kk));
                V_charge{kk,1} = V_(pulse_head:maxlocs(kk));
                T_charge{kk,1} = T_(pulse_head:maxlocs(kk));
                S_charge{kk,1} = S_(pulse_head:maxlocs(kk));
            else 
                t_charge{kk,1} = t_(minlocs(kk-1):maxlocs(kk));
                I_charge{kk,1} = I_(minlocs(kk-1):maxlocs(kk));
                V_charge{kk,1} = V_(minlocs(kk-1):maxlocs(kk));
                T_charge{kk,1} = T_(minlocs(kk-1):maxlocs(kk));
                S_charge{kk,1} = S_(minlocs(kk-1):maxlocs(kk));
            end
        end
    end
else % local minimum comes first
    if length(maxlocs)==length(minlocs) % same number
        for kk = 1:length(minlocs)            
            if kk == 1                            
                for pulse_head = 1:minlocs(1)  % find pulse_head - where the discharge pulse started
                    diffS = S_(pulse_head+1)-S_(pulse_head);
                    if diffS ~= 0            
                    break
                    end 
                end
                t{kk,1} = t_(pulse_head:minlocs(kk));
                I{kk,1} = I_(pulse_head:minlocs(kk));
                V{kk,1} = V_(pulse_head:minlocs(kk));
                T{kk,1} = T_(pulse_head:minlocs(kk));
                S{kk,1} = S_(pulse_head:minlocs(kk));
            else 
                t{kk,1} = t_(maxlocs(kk-1):minlocs(kk));
                I{kk,1} = I_(maxlocs(kk-1):minlocs(kk));
                V{kk,1} = V_(maxlocs(kk-1):minlocs(kk));
                T{kk,1} = T_(maxlocs(kk-1):minlocs(kk));
                S{kk,1} = S_(maxlocs(kk-1):minlocs(kk));
            end           
        end
        t{length(minlocs)+1,1} = t_(maxlocs(length(maxlocs)):length(S_));
        I{length(minlocs)+1,1} = I_(maxlocs(length(maxlocs)):length(S_));
        V{length(minlocs)+1,1} = V_(maxlocs(length(maxlocs)):length(S_));
        T{length(minlocs)+1,1} = T_(maxlocs(length(maxlocs)):length(S_));
        S{length(minlocs)+1,1} = S_(maxlocs(length(maxlocs)):length(S_)); % dischagre segments 
        
        for kk = 1:length(minlocs)                                        % chagre segments
            t_charge{kk,1} = t_(minlocs(kk):maxlocs(kk));
            I_charge{kk,1} = I_(minlocs(kk):maxlocs(kk));
            V_charge{kk,1} = V_(minlocs(kk):maxlocs(kk));
            T_charge{kk,1} = T_(minlocs(kk):maxlocs(kk));
            S_charge{kk,1} = S_(minlocs(kk):maxlocs(kk));
        end
    else % number of local minimums (one) more than that of local maximums
        for kk = 1:length(minlocs)
            if kk==1
                for pulse_head = 1:minlocs(1)  % find pulse_head - where the discharge pulse started
                    diffS = S_(pulse_head+1)-S_(pulse_head);
                    if diffS ~= 0            
                    break
                    end 
                end
                t{kk,1} = t_(pulse_head:minlocs(kk));
                I{kk,1} = I_(pulse_head:minlocs(kk));
                V{kk,1} = V_(pulse_head:minlocs(kk));
                T{kk,1} = T_(pulse_head:minlocs(kk));
                S{kk,1} = S_(pulse_head:minlocs(kk));
            else
                t{kk,1} = t_(maxlocs(kk-1):minlocs(kk));
                I{kk,1} = I_(maxlocs(kk-1):minlocs(kk));
                V{kk,1} = V_(maxlocs(kk-1):minlocs(kk));
                T{kk,1} = T_(maxlocs(kk-1):minlocs(kk));
                S{kk,1} = S_(maxlocs(kk-1):minlocs(kk)); % dischagre segments 
            end
            
            if kk == length(minlocs)                     % chagre segments
                t_charge{kk,1} = t_(minlocs(kk):length(S_));
                I_charge{kk,1} = I_(minlocs(kk):length(S_));
                V_charge{kk,1} = V_(minlocs(kk):length(S_));
                T_charge{kk,1} = T_(minlocs(kk):length(S_));
                S_charge{kk,1} = S_(minlocs(kk):length(S_));
            else
                t_charge{kk,1} = t_(minlocs(kk):maxlocs(kk));
                I_charge{kk,1} = I_(minlocs(kk):maxlocs(kk));
                V_charge{kk,1} = V_(minlocs(kk):maxlocs(kk));
                T_charge{kk,1} = T_(minlocs(kk):maxlocs(kk));
                S_charge{kk,1} = S_(minlocs(kk):maxlocs(kk));
            end
        end 
    end
end
end   
end

