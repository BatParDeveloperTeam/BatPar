function F =  F_easyECMfit(t,I,V,T,S,NRC,R0,flags,param)
    
    N = numel(t);

    OCV1 = param(1);OCV2 = param(2);
    R = param(3:2+NRC); tau = param(3+NRC:2+2*NRC);
    initdata = transpose(reshape(param(3+2*NRC:end),[NRC N])); % param(2+2*NRC:end), NRC*N by 1 vector --> reshape to NRC by N matrix --> transpose to N by NRC matrix

    F = 0;
    
    L2 = @(time,z) sqrt(trapz(time,z.^2));  % anonymous function; trapz - cumulative trapezoidal numerical integration
    
    for kk = 1:N
        
        F = F + L2(t{kk},V{kk} - localsolve(t{kk},I{kk},T{kk},S{kk},OCV1,OCV2,R0,R,tau,initdata(kk,:),flags)); % trapz [(V_experiment - V_model)^2] dt
    end
    
end

