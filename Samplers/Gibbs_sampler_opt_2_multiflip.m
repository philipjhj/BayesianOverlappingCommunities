function [Z,Q, logPs,test_cond_prob_ratio] = Gibbs_sampler_opt_2(A,T,pms,d)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

% clears the persistent variables of the logjoint script
clear logjoint_opt_2_multiflip;

k = 4;

if k>d
    error('the k-flip is bigger than the number of clusters')
end

no_perm = 2^k;

% Initialization of variables for loop
no_obs = size(A,1); no_feats = size(A,2);
z = ones(no_obs,d);
q = ones(no_feats,d);
%NZQ{2} = zeros(no_obs,no_feats);

suff_stats_c = suff_stats_class;

logP = zeros(1,no_perm);
NZQ = cell(1,no_perm);
suff_stats = cell(1,no_perm);
for i = 1:2^k
    suff_stats{i} = suff_stats_c;
    NZQ{i} = zeros(no_obs,no_feats);
end

for l=1:T
    for i = 1:no_feats
        [rand_idx,perms,perm_diff,full_shadow] = k_flip2(z,q(i,:),k);
        dQ = q;
        dQ(i,rand_idx) = perms(1,:);
        prev_q = q;
        %if l == 1 || length(perm_diff{1}) > 0
        [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_2_multiflip(A,z,dQ,pms,suff_stats{1},'q',i,perm_diff{1},full_shadow{1},NZQ{1},prev_q);
        %end
        prev_q = dQ;
        % Permutationerne skal udregnes rigtigt for k > 4 (rekursivt)
        % FIX shadows!
        %
        %
        for flip = 2:no_perm
            dQ(i,rand_idx) = perms(flip,:);
            [logP(flip), NZQ{flip},suff_stats{flip}] = logjoint_opt_2_multiflip(A,z,dQ,pms,suff_stats{flip-1},['q'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_q);
            prev_q = dQ;
        end
        
        P = exp(logP-max(logP));
        idx = find(rand<cumsum(P)/sum(P),1);
        q(i,rand_idx) = perms(idx,:);
        NZQ{1} = NZQ{idx};
        suff_stats{1} = suff_stats{idx};
        prev_q=q;
        
        
        
    end
    
    logPs(l,1) = logP(idx);
    Q{l} = q;
    
    for i = 1:no_obs
        [rand_idx,perms,perm_diff,full_shadow,perm_order] = k_flip2(q,z(i,:),k);
        orig_z = z;
        dZ = z;
        dZ(i,rand_idx) = perms(1,:);
        prev_z = z;
        %if length(perm_diff{1}) > 0
        [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_2_multiflip(A,dZ,q,pms,suff_stats{1},'z',i,perm_diff{1},full_shadow{1},NZQ{1},prev_z);
        %end
        prev_z = dZ;
        
        for flip = 2:no_perm
            dZ(i,rand_idx) = perms(flip,:);
            [logP(flip), NZQ{flip},suff_stats{flip}] = logjoint_opt_2_multiflip(A,dZ,q,pms,suff_stats{flip-1},['z'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_z);
            prev_z = dZ;
        end
        
        P = exp(logP-max(logP));
        idx = find(rand<cumsum(P)/sum(P),1);
        
        z(i,rand_idx) = perms(idx,:);
        NZQ{1} = NZQ{idx};
        suff_stats{1} = suff_stats{idx};
        prev_z = z;
        if 0 < abs(logjoint(A,z,q,pms)-logP(idx))
            logjoint(A,z,q,pms)-logP(idx)
            perm_order
            full_shadow
            rand_idx'
            q'
            nZ = orig_z(i,:);
            nZ(1,rand_idx')=0;
            sum(q(:,[find(nZ)]),2)'>0
            orig_z(i,:)
        end
        
    end
    
    logPs(l,2) = logP(idx);
    Z{l} = z;
end