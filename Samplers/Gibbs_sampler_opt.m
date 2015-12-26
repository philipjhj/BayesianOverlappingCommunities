function [Z,Q] = Gibbs_sampler_opt(A,T,pms,d)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

no_obs = size(A,1); no_feats = size(A,2);
z = ones(no_obs,d);
q = ones(no_feats,d);
NZQ{1} = zeros(no_obs,no_feats);
NZQ{2} = zeros(no_obs,no_feats);
NZQ_active = NZQ{1};

suff_stats_c = suff_stats_class;
suff_stats{1} = suff_stats_c;
suff_stats{2} = suff_stats_c; 


z_prev = z;
q_prev = q;
logP = zeros(1,2);
for l=1:T
    %clear logjoint_opt;
    for i = 1:no_feats
        for j = 1:d
            dQ = q;
            dQ(i,j) = 0;
            [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt(A,z,dQ,pms,suff_stats_c,'q',[i j],NZQ_active,q_prev);
            dQ(i,j) = 1;
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt(A,z,dQ,pms,suff_stats_c,'q',[i j],NZQ_active,q_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            NZQ_active = NZQ{idx};
            suff_stats_c = suff_stats{idx};
            q(i,j) = idx-1;
            q_prev = q;
        end  
    end
    Q{l} = q;
    %clear logjoint_opt;
    for i = 1:no_obs
        for j = 1:d
            dZ = z;
            dZ(i,j) = 0;
            [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt(A,dZ,q,pms,suff_stats_c,'z',[i j],NZQ_active,z_prev);
            dZ(i,j) = 1;
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt(A,dZ,q,pms,suff_stats_c,'z',[i j],NZQ_active,z_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            NZQ_active = NZQ{idx};
            suff_stats_c = suff_stats{idx};
            z(i,j) = idx-1;
            z_prev = z;
        end
    end
    Z{l} = z;
    %fprintf('Sample %d \n',l);
end