function [Z,Q, logPs,pms_all,suff_stat_all,cpu_time,wall_time] = Gibbs_sampler_final(A,T,pms,d,k)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

% clears the persistent variables of the logjoint script
clear LogJoint;

if nargin < 5
    k = 4;
end

if k>d
    error('the k-flip is bigger than the number of clusters')
end

no_perm = 2^k;

% Initialization of variables for loop
no_obs = size(A,1); no_feats = size(A,2);

logPs = nan(T,3);
pms_all = cell(1,T);
suff_stat_all = cell(1,T);

rng('shuffle')
z = randi([0 1],no_obs,d);%ones(no_obs,d);
q = randi([0 1],no_feats,d);%ones(no_feats,d);

suff_stats_c = SuffStatsClass;

logP = zeros(1,no_perm);
NZQ = cell(1,no_perm);
suff_stats = cell(1,no_perm);
for i = 1:2^k
    suff_stats{i} = suff_stats_c;
    NZQ{i} = zeros(no_obs,no_feats);
end


cpu_time = zeros(1,T+1);
wall_time = zeros(1,T);
cpu_time(1,1) = cputime;
for l=1:T
    tic;
    for i = 1:no_feats
        [rand_idx,perms,perm_diff,full_shadow] = kFlip(z,q(i,:),k);
        dQ = q;
        dQ(i,rand_idx) = perms(1,:);
        prev_q = q;
        if l == 1 && i == 1 %|| 0<length(perm_diff{1})
            [logP(1), NZQ{1},suff_stats{1}] = LogJoint(A,z,dQ,pms,suff_stats{1});
            prev_q = dQ;
        elseif ~isempty(perm_diff{1})
            [logP(1), NZQ{1},suff_stats{1}] = LogJoint(A,z,dQ,pms,suff_stats{1},'q',i,perm_diff{1},full_shadow{1},NZQ{1},prev_q);
            prev_q = dQ;
        end
        
        for flip = 2:no_perm
            dQ(i,rand_idx) = perms(flip,:);
            [logP(flip), NZQ{flip},suff_stats{flip}] = LogJoint(A,z,dQ,pms,suff_stats{flip-1},['q'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_q);
            prev_q = dQ;
        end
        P = exp(logP-max(logP));
        idx = find(rand<cumsum(P)/sum(P),1);
        try
            q(i,rand_idx) = perms(idx,:);
        catch
            disp(cumsum(P)/sum(P))
            disp(P)
            disp(idx)
            disp(size(perms,2))
            disp(numel(rand_idx))
            save('error_run_q')
            error('Index problem with q')
        end
        NZQ{1} = NZQ{idx};
        logP(1) = logP(idx);
        suff_stats{1} = suff_stats{idx};
        %prev_q=q;
    end
    
    logPs(l,1) = logP(idx);
    Q{l} = q;
    
    for i = 1:no_obs
        [rand_idx,perms,perm_diff,full_shadow,perm_order] = kFlip(q,z(i,:),k);
        %orig_z = z;
        dZ = z;
        dZ(i,rand_idx) = perms(1,:);
        prev_z = z;
        if ~isempty(perm_diff{1})
            [logP(1), NZQ{1},suff_stats{1}] = LogJoint(A,dZ,q,pms,suff_stats{1},'z',i,perm_diff{1},full_shadow{1},NZQ{1},prev_z);
            prev_z = dZ;
        end
        
        for flip = 2:no_perm
            dZ(i,rand_idx) = perms(flip,:);
            [logP(flip), NZQ{flip},suff_stats{flip}] = LogJoint(A,dZ,q,pms,suff_stats{flip-1},['z'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_z);
            prev_z = dZ;
        end
        
        P = exp(logP-max(logP));
        idx = find(rand<cumsum(P)/sum(P),1);
        logP(1) = logP(idx);
        try
        z(i,rand_idx) = perms(idx,:);
         catch
            disp(cumsum(P)/sum(P))
            disp(P)
            disp(idx)
            disp(size(perms,2))
            disp(numel(rand_idx))
            Error('Index problem with z')
        end
        
        
        NZQ{1} = NZQ{idx};
        suff_stats{1} = suff_stats{idx};
        %prev_z = z;
    end
    logP(1) = logP(idx);
    logPs(l,2) = logP(idx);
    Z{l} = z;
    for hyper_iterations = 1:20
        for i = 1:size(pms,1)
            for j = 1:size(pms,2)
                pm_log=log(pms(i,j));
                pm_log_new=pm_log+0.1*randn;
                pms_new=pms;
                pms_new(i,j)=exp(pm_log_new);
                [logP(2), NZQ{2}, suff_stats{2}] = LogJoint(A,z,q,pms_new,suff_stats{1},'pms',i,j,[],NZQ{1},pms);
                %fprintf(fid,'iteration: %d, idx1: %d, idx2: %d, difference: \t%f\t%f\n',l,i,j,logP(2)-logjoint(A,z,q,pms_new,suff_stats{1}),logP(1)-logjoint(A,z,q,pms));
                %accept_ratio = exp(logP(2))/exp(logP(1));
                P=exp(logP(1:2)-max(logP(1:2)));
                accept_ratio=P(2)/P(1);
                %fprintf('Ratio: %f\n', accept_ratio)
                if rand < accept_ratio
                    %accepted=accepted+1;
                    %fprintf('accepted: %f, change: %f\n',logP(2)/logP(1),pms_new(i,j)-pms(i,j))
                    logP(1)=logP(2);
                    pms=pms_new;
                    suff_stats{1}=suff_stats{2};
                end
            end
        end
        
    end
    %disp(pms)
    pms_all{l} = pms;
    logPs(l,3) = logP(1);
    suff_stat_all{l} = suff_stats{1};
    suff_stat_all{l}.NPap = nan;
    suff_stat_all{l}.NPam = nan;
    suff_stat_all{l}.NPbp = nan;
    suff_stat_all{l}.NPbm = nan;
    
    cpu_time(1,l+1) = cputime;
    wall_time(1,l) = toc;
end
end