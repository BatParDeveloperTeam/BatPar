function [Vhat, t_useful] = fullsolve(t,I,T,S,S_j,OCV,R0,R,tau,flags,NRC) % calculate terminal volatge based on parameterised model and experiment current
% t{kk},I{kk},T{kk},S{kk},         S_j,OCV,R0,R,tau,flags,NRC

% 2nd May 2022 - filter out unavailable SOC
S_j = S_j(~isnan(OCV));
index = ((S >= min(S_j)) & (S <= max(S_j)));  % index of avaible SOC

t = t(index);   %  available timeline
I = I(index);   
T = T(index);   
S = S(index);

OCV = OCV (~isnan(OCV));  % available tables
R0 = R0 (~isnan(R0));
R = reshape(R(~isnan(R)),NRC,[]);
tau = reshape(tau(~isnan(tau)),NRC,[]);

% 29th June 2022 - change from isempty(t) to length(t)<2
if length(t)<2 % no available timeline
    Vhat = nan;
    t_useful = nan;
else

maxdt = 0.001;
tnew = sort(vertcat((t(1):maxdt:t(end))',t(end))); 
tnew = unique(tnew);

I      = interp1(t,I,tnew);
T      = interp1(t,T,tnew);  % I/T/S - column vector
S      = interp1(t,S,tnew);

if (flags.R0temp.tf == 1)                        %
    % reconstruct function for R0                %
    R0 = interpn(flags.R0temp.val,S_j,R0,T,S);   %
else
    R0 = interp1(S_j,R0,S);                      % interploate S_j-R0 to the range of S; R0 - column vector
end

OCV = interp1(S_j,OCV,S);                        % interploate S_j-OCV to the range of S; OCV - column vector



R = nlininterpvec(S_j,R,S);                      % interploate S_j-R to the range of S; R - NRC by numel(S_j) to length(tnew) by NRC --> a hidden transpose!!!

tau = nlininterpvec(S_j,tau,S);                  % interploate S_j-tau to the range of S; tau - NRC by numel(S_j) to length(tnew) by NRC --> a hidden transpose!!!

len = length(tnew);
dt = diff(tnew);

Vi = zeros(len,NRC);

% Integrate up

A = 1./tau;    % tau-length(tnew) by NRC
B = R./tau;    % R-length(tnew) by NRC

% for k = 2:len   
%     Vi(k,:) = ( 1 - 0.5*dt(k-1)*A(k-1,:) ) .* Vi(k-1,:) ...
%                 + (dt(k-1)/2)*B(k-1,:)*(I(k-1) + I(k));
% end

% TYPO FIXED - 17/2/21
for k = 2:len   
    Vi(k,:) = ( 1 - dt(k-1)*A(k-1,:) ) .* Vi(k-1,:) ...
                + (dt(k-1)/2)*B(k-1,:)*(I(k-1) + I(k));
end

Vhat = OCV - R0.*I - sum(Vi,2); % sum over 2nd dim here

Vhat = interp1(tnew,Vhat,t); % could do this faster if we had the indices of the t in the tnew

t_useful = t;
end
end

%{
% Implicit Euler (3RC)
for k = 2:len   
    Vi(k,:) = ( 1 - dt(k-1)*[p(k-1,3),p(k-1,4),p(k-1,5)] ) .* Vi(k-1,:) ...
                + dt(k-1)*[p(k-1,3),p(k-1,4),p(k-1,5)].*[p(k-1,6)*I(k-1),p(k-1,7)*I(k-1),p(k-1,8)*I(k-1)];
end
%}
%{ 
% Midpoint method  (3RC)
for k = 2:len   
    Vi(k,:) = ( 1 - 0.5*dt(k-1)*[p(k-1,3)+p(k,3),p(k-1,4)+p(k,4),p(k-1,5)+p(k,5)] ) .* Vi(k-1,:) ...
                + (dt(k-1)/2)*[p(k-1,3)*p(k-1,6)*I(k-1),p(k-1,4)*p(k-1,7)*I(k-1),p(k-1,5)*p(k-1,8)*I(k-1)] ...
                + (dt(k-1)/2)*[p(k,3)*p(k,6)*I(k),p(k,4)*p(k,7)*I(k),p(k,5)*p(k,8)*I(k)];
end
%}
%{ 
% Midpoint method  (2RC)
for k = 2:len   
    Vi(k,:) = ( 1 - 0.5*dt(k-1)*[p(k-1,3)+p(k,3),p(k-1,4)+p(k,4)] ) .* Vi(k-1,:) ...
                + (dt(k-1)/2)*[p(k-1,3)*p(k-1,5)*I(k-1),p(k-1,4)*p(k-1,6)*I(k-1)] ...
                + (dt(k-1)/2)*[p(k,3)*p(k,5)*I(k),p(k,4)*p(k,6)*I(k)];
end
%}