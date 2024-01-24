% check if tau_i and/or R_i need to be relimited. If yes --> inform the user

tau1_limit = str2num(evalin('base','tau1_limit'));
tau2_limit = str2num(evalin('base','tau2_limit'));
tau3_limit = str2num(evalin('base','tau3_limit'));
Ri_limit = str2num(evalin('base','Ri_limit'));

if NRC == 1
   if isempty(tau1_limit) == 1
        if tau_1_reco >= 1000 * 0.95 
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_1_reco >= tau1_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end 

   RR1 = prctile (R(1,:),50);   % median of R1
   
   if isempty(Ri_limit) == 1
        if RR1 >= 0.3 * 0.95 
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if RR1 >= Ri_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end
end

if NRC == 2
   if isempty(tau1_limit) == 1
        if tau_1_reco >= 50 * 0.95 
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_1_reco >= tau1_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end 
   
   if isempty(tau2_limit) == 1
        if tau_2_reco >= 1500 * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau2 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_2_reco >= tau2_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau2 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end 

   RR1 = prctile (R(1,:),50);   % median of R1
   RR2 = prctile (R(2,:),50);   % median of R2

   if isempty(Ri_limit) == 1
        if (RR1 >= 0.3 * 0.95) || (RR2 >= 0.3 * 0.95) 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if (RR1 >= Ri_limit * 0.95) || (RR2 >= Ri_limit * 0.95)
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end
end

if NRC == 3
   if isempty(tau1_limit) == 1
        if tau_1_reco >= 10 * 0.95 
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_1_reco >= tau1_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau1 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end 
   
   if isempty(tau2_limit) == 1
        if tau_2_reco >= 100 * 0.95 
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau2 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_2_reco >= tau2_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau2 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end 

   if isempty(tau3_limit) == 1
        if tau_3_reco >= 1000 * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau3 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if tau_3_reco >= tau3_limit * 0.95 
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for tau3 and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end

   RR1 = prctile (R(1,:),50);   % median of R1
   RR2 = prctile (R(2,:),50);   % median of R2
   RR3 = prctile (R(3,:),50);   % median of R3 

   if isempty(Ri_limit) == 1
        if (RR1 >= 0.3 * 0.95) || (RR2 >= 0.3 * 0.95) || (RR3 >= 0.3 * 0.95 )
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   else 
        if (RR1 >= Ri_limit * 0.95) || (RR2 >= Ri_limit * 0.95) || (RR3 >= Ri_limit * 0.95)
           disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
            disp('     Please increase the upper limit for Ri and rerun the codes')
            disp('     ======CAUTION======CAUTION======CAUTION======CAUTION======')
        end
   end
end
