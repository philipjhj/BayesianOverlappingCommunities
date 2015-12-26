function [Z,Q, logPs] = Gibbs_sampler_opt_2(A,T,pms,d)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

clear logjoint_opt_num_int;

no_obs = size(A,1); no_feats = size(A,2);
z = ones(no_obs,d);
q = ones(no_feats,d);
NZQ{1} = zeros(no_obs,no_feats);
NZQ{2} = zeros(no_obs,no_feats);

suff_stats_c = suff_stats_class;
suff_stats{1} = suff_stats_c;
suff_stats{2} = suff_stats_c; 

z_prev = z;
q_prev = q;

logP = zeros(1,2);
for l=1:T
    [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_num_int(A,z,q,pms,suff_stats{1},'q',[1 1],NZQ{1},q);
    for i = 1:no_feats
        for j = 1:d
            dQ = q;
            dQ(i,j) = abs(dQ(i,j)-1);
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt_num_int(A,z,dQ,pms,suff_stats{1},'q',[i j],NZQ{1},q_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            logP(1) = logP(idx);
            NZQ{1} = NZQ{idx};
            suff_stats{1} = suff_stats{idx};
            q(i,j) = dQ(i,j) == (idx-1);
            q_prev = q;
        end
    end
    logPs(l,1) = logP(idx);
    Q{l} = q;
    [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_num_int(A,z,q,pms,suff_stats{1},'z',[1 1],NZQ{1},z);
    for i = 1:no_obs        
        for j = 1:d
            dZ = z;
            dZ(i,j) = abs(dZ(i,j)-1);
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt_num_int(A,dZ,q,pms,suff_stats{1},'z',[i j],NZQ{1},z_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            logP(1) = logP(idx);
            NZQ{1} = NZQ{idx};
            suff_stats{1} = suff_stats{idx};
            z(i,j) = dZ(i,j) == (idx-1);
            z_prev = z;
        end
    end
    logPs(l,2) = logP(idx);
    Z{l} = z;
    %fprintf('Sample %d \n',l);
end