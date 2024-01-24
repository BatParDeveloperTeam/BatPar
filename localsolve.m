function Vhat = localsolve(t,I,T,S,OCV1,OCV2,R0,R,tau,initdata,flags)

% Mapping the data from old timeline to new, fixed-step timeline (tnew)
% The resaon for doing this is that experimentalists sometimes use variable time steps, in which case the data is not evenly logged with respect to time
% However, the ODE solver runs on the basis of time logged, so uneven data is likey to casue local underfit or local overfit over time

maxdt = 0.05;  
tnew = sort(vertcat((t(1):maxdt:t(end))',t(end))); % subdivide timeline; vertcat - splice two vectors; sort - ascending order
tnew = unique(tnew);                               % unique tnew in case t(end) can be divisible by maxdt, in which case there would be two 't(end)' in tnew 

I      = interp1(t,I,tnew);         
T      = interp1(t,T,tnew);        % interpolate new data on tnew
S      = interp1(t,S,tnew);

if (flags.R0temp.tf == 1)  % 
    R0 = R0(T);            %
end                        %

len = length(tnew);  
dt = diff(tnew);          

Vi = zeros(len,length(initdata));    % Voltage of each RC pair; initdata - N by NRC vector, length(initdata) = NRC
                                     % Vi - len by NRC matrix 

Vi(1,:) = initdata + 1E-10*sort(rand(size(initdata)),'descend'); % random noise to split tau_1 and tau_2, R1 and R2. Otherwise, the optimiser would probably give tau_1=tau_2 and R1=R2 
% note, initdata is auto-transposed by this assignment

% Integrate up; See Ruben's slides --> 3RC ECM

A = transpose(1./tau);   % tau - NRC by 1 vector
B = transpose(R./tau);   % column vector transposed to row vector

% Calculate Vi on tnew
% Vi_(k) = Vi_(k-1) + delta[Vi]_(k); delta[Vi]_(k) ~= -dt_(k-1)* [A*Vi_(k-1) - B*Ii_(k-1)]

for k = 2:len   
    Vi(k,:) = ( 1 - dt(k-1)*A ) .* Vi(k-1,:) ...
                + (dt(k-1)/2)*B*(I(k-1) + I(k)); % (I(k-1) + I(k))/2 - midpoint
end

% for k = 2:len   
%     Vi(k,:) = ( 1 - 0.5*dt(k-1)*A ) .* Vi(k-1,:) ...
%                 + (dt(k-1)/2)*B*(I(k-1) + I(k));
% end

Vhat = OCV1*S+OCV2- R0.*I - sum(Vi,2);  % OCV & R0 - single number; 
                                           % S & I - column vector (the vector holding data at around S_j); 
                                           % sum(Vi,2) - sum up the voltage of each RC pair; '2' means sum over the 2 dimsension
                                           % Vhat - column vector
Vhat = interp1(tnew,Vhat,t);   % interpolate back onto t

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