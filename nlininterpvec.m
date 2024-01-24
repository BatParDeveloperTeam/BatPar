function outvec = nlininterpvec(points,theta,evaluateat)   % (S_j,R,S) or (S_j,tau,S)    
% points - row vector  theta - NRC by numel(S_j) evaluateat - column vector

N = size(theta,1);     % N= NRC; 
outvec = zeros(size(evaluateat,1),N);    % outvec - length(S) by NRC; length(S)=length(tnew)

for ii = 1:N
    outvec(:,ii) = interp1(points,theta(ii,:),evaluateat); % outvec(:,ii) - column vector; here is a hidden transpose from!!!
end

end

