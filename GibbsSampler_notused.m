function myGibbsChain = GibbsSampler(myGibbsChain,testError,ActiveGibbs,NotConstrained,activeHyper,hyperindex)
% Gibbs sampler
% A is data
% T number of samples
% pms parameters
% d number of concepts

% clears the persistent variables of the logjoint script
clear LogJoint;

if nargin<2
    testError = true;
end
if nargin<6
    hyperindex = [];
end

if myGibbsChain.k>myGibbsChain.d
    error('the k-flip is bigger than the number of clusters')
end
no_perm = 2^myGibbsChain.k;
% Initialization of variables for loop
nObs = size(myGibbsChain.A,1); nFeats = size(myGibbsChain.A,2);

%logPs = nan(T,3);
%pms_all = cell(1,T);
%suff_stat_all = cell(1,T);
%rng(1)
%rng('shuffle')
if isempty(myGibbsChain.zSamples{1})
    
    z = randi([0 1],nObs,myGibbsChain.d);%ones(no_obs,d);
    q = randi([0 1],nFeats,myGibbsChain.d);%ones(no_feats,d);
else
    lastIter = find(~cellfun('isempty',myGibbsChain.zSamples),1,'last');
    z = myGibbsChain.zSamples{lastIter};
    q = myGibbsChain.qSamples{lastIter};
end


suff_stats_c = SuffStatsClass;
logP = zeros(1,no_perm);
NZQ = cell(1,no_perm);
suff_stats = cell(1,no_perm);
for i = 1:2^myGibbsChain.k
    suff_stats{i} = suff_stats_c;
    NZQ{i} = zeros(nObs,nFeats);
end

%Q = cell(1,T);
%Z = cell(1,T);

%cpu_time = zeros(1,T+1);
%wall_time = zeros(1,T);
myGibbsChain.cpuTime(1,1) = cputime;
for l=myGibbsChain.Tstart:myGibbsChain.T
    tic;
    if ActiveGibbs
        for i = 1:nFeats
            [rand_idx,perms,perm_diff,full_shadow] = kFlip(z,q(i,:),myGibbsChain.k);
            dQ = q;
            dQ(i,rand_idx) = perms(1,:);
            prev_q = q;
            if l == 1 && i == 1 %|| 0<length(perm_diff{1})
                [logP(1), NZQ{1},suff_stats{1}] = LogJoint(NotConstrained,myGibbsChain.A,z,dQ,myGibbsChain.pms,suff_stats{1});
                prev_q = dQ;
            elseif ~isempty(perm_diff{1})
                [logP(1), NZQ{1},suff_stats{1}] = LogJoint(NotConstrained,myGibbsChain.A,z,dQ,myGibbsChain.pms,suff_stats{1},'q',i,perm_diff{1},full_shadow{1},NZQ{1},prev_q);
                prev_q = dQ;
            end
            
            for flip = 2:no_perm
                dQ(i,rand_idx) = perms(flip,:);
                [logP(flip), NZQ{flip},suff_stats{flip}] = LogJoint(NotConstrained,myGibbsChain.A,z,dQ,myGibbsChain.pms,suff_stats{flip-1},['q'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_q);
                if isnan(logP(flip))
                    disp(cumsum(P)/sum(P))
                    disp(P)
                    disp(idx)
                    disp(size(perms,2))
                    disp(numel(rand_idx))
                    save('error_run_q')
                    display('Index problem with q')
                    return
                end
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
        
        myGibbsChain.logPs(1,l) = logP(idx);
        myGibbsChain.qSamples{l} = q;
        
        if testError
            mode = 'q';
            test_suff=calc_suff_stats(myGibbsChain.A,z,q,myGibbsChain.pms);
            status=test_suff_stats(suff_stats{1},test_suff,mode,myGibbsChain.pms);
            if status == -1
                myGibbsChain.Tstart=l+1;
                save(strcat('Output/ErrorRuns/error_',num2str(now),'.mat'))
                display('suff_stat is wrong, data saved to output folder')
                display('Ending Gibbs Sampler')
                return
            end
        end
        
        for i = 1:nObs
            [rand_idx,perms,perm_diff,full_shadow,perm_order] = kFlip(q,z(i,:),myGibbsChain.k);
            %orig_z = z;
            dZ = z;
            dZ(i,rand_idx) = perms(1,:);
            prev_z = z;
            if ~isempty(perm_diff{1})
                [logP(1), NZQ{1},suff_stats{1}] = LogJoint(NotConstrained,myGibbsChain.A,dZ,q,myGibbsChain.pms,suff_stats{1},'z',i,perm_diff{1},full_shadow{1},NZQ{1},prev_z);
                prev_z = dZ;
            end
            
            for flip = 2:no_perm
                dZ(i,rand_idx) = perms(flip,:);
                [logP(flip), NZQ{flip},suff_stats{flip}] = LogJoint(NotConstrained,myGibbsChain.A,dZ,q,myGibbsChain.pms,suff_stats{flip-1},['z'],i,perm_diff{flip},full_shadow{flip},NZQ{flip-1},prev_z);
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
        myGibbsChain.logPs(2,l) = logP(idx);
        myGibbsChain.zSamples{l} = z;
        
        if testError
            mode = 'z';
            test_suff=calc_suff_stats(myGibbsChain.A,z,q,myGibbsChain.pms);
            status=test_suff_stats(suff_stats{1},test_suff,mode,myGibbsChain.pms);
            if status == -1
                myGibbsChain.Tstart=l+1;
                save(strcat('Output/ErrorRuns/error_',num2str(now),'.mat'))
                display('suff_stat is wrong, data saved to output folder')
                display('Ending Gibbs Sampler')
                return
            end
        end
    end
    if activeHyper
        firstDone=false;
        for hyper_iterations = 1:20
            for i = 1:size(myGibbsChain.pms,1)
                for j = 1:size(myGibbsChain.pms,2)
                    if ~isempty(hyperindex)
                        if any(i == hyperindex(1,:)) || any(j == hyperindex(2,:));
                            continue
                        end
                    end
                    pm_log=log(myGibbsChain.pms(i,j));
                    pm_log_new=pm_log+0.1*randn;
                    pms_new=myGibbsChain.pms;
                    pms_new(i,j)=exp(pm_log_new);
                    if l == myGibbsChain.Tstart && ActiveGibbs==false && ~firstDone %|| 0<length(perm_diff{1})
                        firstDone = true;
                        [logP(1), NZQ{1},suff_stats{1}] = LogJoint(NotConstrained,myGibbsChain.A,z,q,myGibbsChain.pms,suff_stats{1});
                    end
                    [logP(2), NZQ{2}, suff_stats{2}] = LogJoint(NotConstrained,myGibbsChain.A,z,q,pms_new,suff_stats{1},'pms',i,j,[],NZQ{1},myGibbsChain.pms);
                    %fprintf(fid,'iteration: %d, idx1: %d, idx2: %d, difference: \t%f\t%f\n',l,i,j,logP(2)-logjoint(A,z,q,pms_new,suff_stats{1}),logP(1)-logjoint(A,z,q,pms));
                    %accept_ratio = exp(logP(2))/exp(logP(1));
                    P=exp(logP(1:2)-max(logP(1:2)));
                    %display(P);
                    accept_ratio=P(2)/P(1)*pms_new(i,j)/myGibbsChain.pms(i,j);
                    %fprintf('Ratio: %f\n', accept_ratio)
                    if rand < accept_ratio
                        %accepted=accepted+1;
                        %fprintf('accepted: %f, change: %f\n',logP(2)/logP(1),pms_new(i,j)-pms(i,j))
                        logP(1)=logP(2);
                        myGibbsChain.pms=pms_new;
                        suff_stats{1}=suff_stats{2};
                    end
                end
            end
            
        end
        %disp(pms)
        myGibbsChain.pmsSamples{l} = myGibbsChain.pms;
        myGibbsChain.logPs(3,l) = logP(1);
        myGibbsChain.SuffStatAll{l} = suff_stats{1};
        myGibbsChain.SuffStatAll{l}.NPap = nan;
        myGibbsChain.SuffStatAll{l}.NPam = nan;
        myGibbsChain.SuffStatAll{l}.NPbp = nan;
        myGibbsChain.SuffStatAll{l}.NPbm = nan;
        
        if testError
            mode = 'pms';
            test_suff=calc_suff_stats(myGibbsChain.A,z,q,myGibbsChain.pms);
            status=test_suff_stats(suff_stats{1},test_suff,mode,myGibbsChain.pms);
            if status == -1
                myGibbsChain.Tstart=l+1;
                save(strcat('Output/ErrorRuns/error_',num2str(now),'.mat'))
                display('suff_stat is wrong, data saved to output folder')
                display('Ending Gibbs Sampler')
                return
            end
        end
    end
    myGibbsChain.cpuTime(1,l+1) = cputime;
    myGibbsChain.wallTime(1,l) = toc;
    fprintf('Sample: %d/%d - Time: %f \n',l,myGibbsChain.T,sum(myGibbsChain.wallTime(1,1:l)));
end


end

function status = test_suff_stats(suff_stat_step,suff_stat_true,mode,pms)

fields = properties(suff_stat_true);
for i=1:numel(fields)
    if numel(suff_stat_true.(fields{i})) > 0
        %disp(fields{i})
        diff=sum(suff_stat_step.(fields{i})-suff_stat_true.(fields{i}));
        
        if abs(diff)>1e-6
            disp('Error with suff_stat')
            disp(fields{i})
            disp('mode')
            disp(mode)
            disp('diff')
            disp(suff_stat_step.(fields{i})-suff_stat_true.(fields{i}))
            disp('hyperparameters')
            disp(pms)
            %             if mode(1) == 'p'
            %                 disp(pms_old)
            %             end
            status = -1;
            return
        end
    end
end
status=1;
end

function suff_stats = calc_suff_stats(A,Z,Q,pms)
no_obs = size(A,1);
no_feats = size(A,2);
no_concepts = size(Z,2);

suff_stats = SuffStatsClass;
% #Concepts shared for observation i and feature j
NZQ = Z*Q';

% #Elements that share concepts and are linked
suff_stats.NPap = A.*(0<NZQ);
% #Elements that share concepts and aren't linked
suff_stats.NPam = (1-A).*(0<NZQ);
% #Elements that does not share concepts and are linked
suff_stats.NPbp = A.*(1-(0<NZQ));
% #Elements that do not share concepts and aren't linked
suff_stats.NPbm = (1-A).*(1-(0<NZQ));

% prints to check for handling of missing data
%[sum(suff_stats.NPap);
%sum(suff_stats.NPam);
%sum(suff_stats.NPbp);
%sum(suff_stats.NPbm)]'

% Handle missing data (NaNs)
suff_stats.NPap(isnan(suff_stats.NPap)) = 0;
suff_stats.NPam(isnan(suff_stats.NPam)) = 0;
suff_stats.NPbp(isnan(suff_stats.NPbp)) = 0;
suff_stats.NPbm(isnan(suff_stats.NPbm)) = 0;

suff_stats.NZp = sum(Z)'+ pms(1,1);
suff_stats.NZm = no_obs-sum(Z)'+ pms(1,2);

suff_stats.NQp = sum(Q)'+pms(2,1);
suff_stats.NQm = no_feats-sum(Q)'+pms(2,2);

% prints to check for handling of missing data
%[sum(suff_stats.NPap);
%sum(suff_stats.NPam);
%sum(suff_stats.NPbp);
%sum(suff_stats.NPbm)]'

% Adding the parameters to get
% final input for beta functions
suff_stats.NPaps = sum(suff_stats.NPap)'+pms(3,1);
suff_stats.NPams = sum(suff_stats.NPam)'+pms(3,2);
suff_stats.NPbps = sum(suff_stats.NPbp)'+pms(4,1);
suff_stats.NPbms = sum(suff_stats.NPbm)'+pms(4,2);

%status=logjointNaive(A,Z,Q,pms,suff_stats);
%status = -1 means problem with sufficient stats, 1 for ok
end