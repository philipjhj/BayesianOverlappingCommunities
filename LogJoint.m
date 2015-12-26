function [logprob,NZQ,suff_stats] = logjoint_final(A,Z,Q,pms,suff_stats,mode,iter,perm_diff,shadow,NZQ,M_prev)
% Given a matrix A with its corresponding observation X concept matrix
% Z and feature X concept matrix Q, this function calculates
% the probability of this combination given the parameters.

if nargin < 6
    mode = 0;
    idx = 0;
    NZQ = zeros(size(A,1),size(A,2));
    iter=0;
    perm_diff=[];
    shadow=[];
    M_prev=[];
end

persistent no_obs
persistent no_feats
persistent no_concepts

if isempty(no_obs)
    no_obs = size(A,1);
    no_feats = size(A,2);
    no_concepts = size(Z,2);
end


if isempty(suff_stats.NPap)
    
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
    
    
    %suff_stats.betasuff_AZQ = betaln(suff_stats.NPaps,suff_stats.NPams)+betaln(suff_stats.NPbps,suff_stats.NPbms);
    suff_stats.betasuff_AZQ = betaintegral_noverbose(suff_stats.NPaps,suff_stats.NPams,suff_stats.NPbps,suff_stats.NPbms,50000);
    suff_stats.betahyper_probs = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));
    
    suff_stats.betasuff_sums = betaln(suff_stats.NZp,suff_stats.NZm)+betaln(suff_stats.NQp,suff_stats.NQm);
    suff_stats.betahyper_ZQ = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));
    
    
elseif mode(1) == 'z'
    %change = str2num(mode(2));
    i = iter;
    d = perm_diff;
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Z(i,d) is changed.
    %      NZQ(i,:)    +  Z(i,d)*Q(:,d)' - Z_prev(i,d)*Q(:,d)'
    % (previous value) + (change compared to previous value)
    
    Z_diff = Z(i,d)-M_prev(i,d);
    if Z_diff ~= 0
        NZQ(i,:) = NZQ(i,:)+sum(repmat(Z_diff,size(Q,1),1)'.*Q(:,d)',1);
    end
    
    % Find difference from previous state
    suff_stats.NZp(d) = suff_stats.NZp(d)+Z_diff';
    suff_stats.NZm(d) = suff_stats.NZm(d)-Z_diff';
    
    
    if length(shadow) > 0
        %shadow = 1:size(suff_stats.NPap(i,:),2);
        NPap_old = suff_stats.NPap(i,shadow);
        NPam_old = suff_stats.NPam(i,shadow);
        NPbp_old = suff_stats.NPbp(i,shadow);
        NPbm_old = suff_stats.NPbm(i,shadow);
        
        PZQ =(0<(NZQ(i,shadow)));
        AZQ = A(i,shadow).*PZQ;
        
        % #Elements that share concepts and are linked
        suff_stats.NPap(i,shadow) = AZQ;
        % #Elements that share concepts and aren't linked
        suff_stats.NPam(i,shadow) = PZQ-AZQ;
        % #Elements that does not share concepts and are linked
        suff_stats.NPbp(i,shadow) = A(i,shadow)-AZQ;
        % #Elements that do not share concepts and aren't linked
        suff_stats.NPbm(i,shadow) = AZQ+1-PZQ-A(i,shadow); %(1-A(i,:)).*(1-(0<(NZQ(i,:))));
        
        % Handle missing data (NaNs)
        suff_stats.NPap(isnan(suff_stats.NPap)) = 0;
        suff_stats.NPam(isnan(suff_stats.NPam)) = 0;
        suff_stats.NPbp(isnan(suff_stats.NPbp)) = 0;
        suff_stats.NPbm(isnan(suff_stats.NPbm)) = 0;
        
        
        suff_stats.NPaps(shadow) = suff_stats.NPaps(shadow)+suff_stats.NPap(i,shadow)'-NPap_old';
        suff_stats.NPams(shadow) = suff_stats.NPams(shadow)+suff_stats.NPam(i,shadow)'-NPam_old';
        suff_stats.NPbps(shadow) = suff_stats.NPbps(shadow)+suff_stats.NPbp(i,shadow)'-NPbp_old';
        suff_stats.NPbms(shadow) = suff_stats.NPbms(shadow)+suff_stats.NPbm(i,shadow)'-NPbm_old';
        
        suff_stats.betasuff_AZQ(shadow) = betaintegral_noverbose(suff_stats.NPaps(shadow),suff_stats.NPams(shadow),suff_stats.NPbps(shadow),suff_stats.NPbms(shadow),50000);
        %suff_stats.betasuff_AZQ(shadow) = betaln(suff_stats.NPaps(shadow),suff_stats.NPams(shadow))+betaln(suff_stats.NPbps(shadow),suff_stats.NPbms(shadow));
    end
    %suff_stats.betahyper_probs = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));
    
    suff_stats.betasuff_sums(d) = betaln(suff_stats.NZp(d),suff_stats.NZm(d))+betaln(suff_stats.NQp(d),suff_stats.NQm(d));
    %suff_stats.betahyper_ZQ = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));
    
    
elseif mode(1) == 'q'
    j = iter;
    d = perm_diff;
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Q(j,d) is changed.
    %      NZQ(:,j)    +  Z(:,d)*Q(j,d)' - Z(:,d)*Q_prev(j,d)'
    % (previous value) + (change compared to previous value)
    
    Q_diff=Q(j,d)-M_prev(j,d);
    if Q_diff ~= 0
        NZQ(:,j) = NZQ(:,j)+sum(Z(:,d).*repmat(Q_diff,size(Z,1),1),2);
    end
    
    
    % Find difference from previous state
    suff_stats.NQp(d) = suff_stats.NQp(d)+Q_diff';
    suff_stats.NQm(d) = suff_stats.NQm(d)-Q_diff';
    
    if length(shadow) > 0
        %    shadow = 1:size(suff_stats.NPap(:,j),1);
        
        NPap_old = suff_stats.NPap(shadow,j);
        NPam_old = suff_stats.NPam(shadow,j);
        NPbp_old = suff_stats.NPbp(shadow,j);
        NPbm_old = suff_stats.NPbm(shadow,j);
        
        PZQ =(0<(NZQ(shadow,j)));
        AZQ = A(shadow,j).*PZQ;
        
        % #Elements that share concepts and are linked
        suff_stats.NPap(shadow,j) = AZQ;
        % #Elements that share concepts and aren't linked
        suff_stats.NPam(shadow,j) = PZQ-AZQ;
        % #Elements that does not share concepts and are linked
        suff_stats.NPbp(shadow,j) = A(shadow,j)-AZQ;
        % #Elements that do not share concepts and aren't linked
        suff_stats.NPbm(shadow,j) = AZQ+1-PZQ-A(shadow,j); %(1-A(i,:)).*(1-(0<(NZQ(i,:))));
        
        % Handle missing data (NaNs)
        suff_stats.NPap(isnan(suff_stats.NPap)) = 0;
        suff_stats.NPam(isnan(suff_stats.NPam)) = 0;
        suff_stats.NPbp(isnan(suff_stats.NPbp)) = 0;
        suff_stats.NPbm(isnan(suff_stats.NPbm)) = 0;
        
        
        suff_stats.NPaps(j) = suff_stats.NPaps(j)+sum(suff_stats.NPap(shadow,j)-NPap_old);
        suff_stats.NPams(j) = suff_stats.NPams(j)+sum(suff_stats.NPam(shadow,j)-NPam_old);
        suff_stats.NPbps(j) = suff_stats.NPbps(j)+sum(suff_stats.NPbp(shadow,j)-NPbp_old);
        suff_stats.NPbms(j) = suff_stats.NPbms(j)+sum(suff_stats.NPbm(shadow,j)-NPbm_old);
        
        
        suff_stats.betasuff_AZQ(j) = betaintegral_noverbose(suff_stats.NPaps(j),suff_stats.NPams(j),suff_stats.NPbps(j),suff_stats.NPbms(j),50000);
        %suff_stats.betasuff_AZQ(j) = betaln(suff_stats.NPaps(j),suff_stats.NPams(j))+betaln(suff_stats.NPbps(j),suff_stats.NPbms(j));
    end
    
    try
    suff_stats.betasuff_sums(d) = betaln(suff_stats.NZp(d),suff_stats.NZm(d))+betaln(suff_stats.NQp(d),suff_stats.NQm(d));
    catch
        disp(suff_stats.NZp(d))
        disp(suff_stats.NZm(d))
        disp(suff_stats.NQp(d))
        disp(suff_stats.NQm(d))
        disp(sum(Z(:,d).*repmat(Q_diff,size(Z,1),1),2))
    end
    
    
    
elseif mode(1) == 'p'
    pms_old=M_prev;
    i = iter;
    j = perm_diff;
    suff_stat_change = pms(i,j)-pms_old(i,j);
    if i == 2 && j == 1
        suff_stats.NQp = suff_stats.NQp+suff_stat_change;
    elseif i == 2 && j == 2
        suff_stats.NQm = suff_stats.NQm+suff_stat_change;
    elseif i == 1 && j == 1
        suff_stats.NZp = suff_stats.NZp+suff_stat_change;
    elseif i == 1 && j == 2
        suff_stats.NZm = suff_stats.NZm+suff_stat_change;
    elseif i == 3 && j == 1
        suff_stats.NPaps = suff_stats.NPaps+suff_stat_change;
    elseif i == 3 && j == 2
        suff_stats.NPams = suff_stats.NPams+suff_stat_change;
    elseif i == 4 && j == 1
        suff_stats.NPbps = suff_stats.NPbps+suff_stat_change;
    elseif i == 4 && j == 2
        suff_stats.NPbms = suff_stats.NPbms + suff_stat_change;
    end
    
    
    
    
    if any(i == [1 2])
        try
            suff_stats.betasuff_sums = betaln(suff_stats.NZp,suff_stats.NZm)+betaln(suff_stats.NQp,suff_stats.NQm);
        catch
            disp('NZp')
            disp(suff_stats.NZp)
            disp('NZm')
            disp(suff_stats.NZm)
            disp('NQp')
            disp(suff_stats.NQp)
            disp('NQm')
            disp(suff_stats.NQm)
            disp(pms_old(i,j))
            disp(pms(i,j))
            save('error_run_beta')
            error('non negative input to betaln()')
        end
        
        suff_stats.betahyper_ZQ = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));
    else
        suff_stats.betahyper_probs = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));
        %suff_stats.betasuff_AZQ = betaln(suff_stats.NPaps,suff_stats.NPams)+betaln(suff_stats.NPbps,suff_stats.NPbms);
        suff_stats.betasuff_AZQ = betaintegral_noverbose(suff_stats.NPaps,suff_stats.NPams,suff_stats.NPbps,suff_stats.NPbms,50000);
    end
end
% 
% disp(suff_stats.NPaps)
% disp(suff_stats.NPams)
% disp(suff_stats.NPbps)
% disp(suff_stats.NPbms)


test_suff=calc_suff_stats(A,Z,Q,pms);
test_suff_stats(suff_stats,test_suff,mode,pms,M_prev)
%disp(test_suff)

%disp(betaintegral_noverbose(suff_stats.NPaps(5),suff_stats.NPams(5),suff_stats.NPbps(5)+1,suff_stats.NPbms(5),50000))

%disp(suff_stats.betahyper_ZQ)
%disp(suff_stats.betahyper_probs)
%disp(suff_stats.betasuff_AZQ)
%disp(suff_stats.betasuff_sums)
% Calculating the probability
logprob=sum(suff_stats.betasuff_AZQ)-no_feats*suff_stats.betahyper_probs+sum(suff_stats.betasuff_sums)-no_concepts*suff_stats.betahyper_ZQ;

end


function test_suff_stats(suff_stat_step,suff_stat_true,mode,pms,pms_old)

fields = properties(suff_stat_true);
for i=1:numel(fields)
    if numel(suff_stat_true.(fields{i})) > 0
        %disp(fields{i})
        diff=sum(sum(suff_stat_step.(fields{i})-suff_stat_true.(fields{i})));
        if abs(diff)>1e-10
            disp('Error')
            disp(fields{i})
            disp(mode)
            disp(suff_stat_step.(fields{i})-suff_stat_true.(fields{i}))
            disp(pms)
            if mode(1) == 'p'
                disp(pms_old)
            end
            error('suff_stat is wrong')
        end
    end
end




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

end
