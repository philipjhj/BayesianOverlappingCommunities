function [Z,Q, logPs,pms_all] = Gibbs_sampler_opt_2_hyper(A,T,pms,d)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

clear logjoint_opt_2_hyper;

no_obs = size(A,1); no_feats = size(A,2);
z = ones(no_obs,d);
q = ones(no_feats,d);
NZQ{1} = zeros(no_obs,no_feats);
NZQ{2} = zeros(no_obs,no_feats);

suff_stats_c = suff_stats_class;
suff_stats{1} = suff_stats_c;
suff_stats{2} = suff_stats_c; 

pms_all = cell(T,1);
% accepted=0;

z_prev = z;
q_prev = q;
%delete('logpdiff.txt');
%fid = fopen('logpdiff.txt','a');

logP = zeros(1,2);
logPs = zeros(T,3);
for l=1:T
    [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_2_hyper(A,z,q,pms,suff_stats{1},'q',[1 1],NZQ{1},q);
    for i = 1:no_feats
        for j = 1:d
            dQ = q;
            dQ(i,j) = abs(dQ(i,j)-1);
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt_2_hyper(A,z,dQ,pms,suff_stats{1},'q',[i j],NZQ{1},q_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            logP(1) = logP(idx);
            NZQ{1} = NZQ{idx};
            suff_stats{1} = suff_stats{idx};
            q(i,j) = dQ(i,j) == (idx-1);
            q_prev = q;
        end
    end
    logPs(l,1) = logP(1);
    Q{l} = q;
    [logP(1), NZQ{1},suff_stats{1}] = logjoint_opt_2_hyper(A,z,q,pms,suff_stats{1},'z',[1 1],NZQ{1},z);
    for i = 1:no_obs        
        for j = 1:d
            dZ = z;
            dZ(i,j) = abs(dZ(i,j)-1);
            [logP(2), NZQ{2},suff_stats{2}] = logjoint_opt_2_hyper(A,dZ,q,pms,suff_stats{1},'z',[i j],NZQ{1},z_prev);
            P = exp(logP-max(logP));
            idx = find(rand<cumsum(P)/sum(P),1);
            logP(1) = logP(idx);
            NZQ{1} = NZQ{idx};
            suff_stats{1} = suff_stats{idx};
            z(i,j) = dZ(i,j) == (idx-1);
            z_prev = z;
        end
    end
    logPs(l,2) = logP(1);
    Z{l} = z;
    
    for i = 1:size(pms,1)
        for j = 1:size(pms,2)
            pm_log=log(pms(i,j));
            pm_log_new=pm_log+0.1*randn;
            pms_new=pms;
            pms_new(i,j)=exp(pm_log_new);
            [logP(2), NZQ{2}, suff_stats{2}] = logjoint_opt_2_hyper(A,z,q,pms_new,suff_stats{1},'pms',[i j],NZQ{1},pms);
            %fprintf(fid,'iteration: %d, idx1: %d, idx2: %d, difference: \t%f\t%f\n',l,i,j,logP(2)-logjoint(A,z,q,pms_new,suff_stats{1}),logP(1)-logjoint(A,z,q,pms));
            P=exp(logP-max(logP));
            %accept_ratio = exp(logP(2))/exp(logP(1));
            accept_ratio = P(2)/P(1);
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
    pms_all{l} = pms;
    logPs(l,3) = logP(1);
%     fprintf('Sample %d \n',l);
end
% accepted
%fclose(fid);
end