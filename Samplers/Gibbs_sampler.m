function [Z,Q] = Gibbs_sampler(A,T,pms,d)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

no_obs = size(A,1); no_feats = size(A,2);
z = ones(no_obs,d);
q = ones(no_feats,d);

for l=1:T
    for i = 1:no_feats
        for j = 1:d
            logP = zeros(1,2);
            dq = q;
            dq(i,j) = 0;
            logP(1) = logjoint(A,z,dq,pms);
            dq(i,j) = 1;
            logP(2) = logjoint(A,z,dq,pms);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            q(i,j) = idx-1;
        end  
    end
    Q{l} = q;

    for i = 1:no_obs
        for j = 1:d
            logP = zeros(1,2);
            dz = z;
            dz(i,j) = 0;
            logP(1) = logjoint(A,dz,q,pms);
            dz(i,j) = 1;
            logP(2) = logjoint(A,dz,q,pms);
            P = exp(logP-max(logP));
            
            idx = find(rand<cumsum(P)/sum(P),1);
            z(i,j) = idx-1;
        end
    end
    Z{l} = z;
end