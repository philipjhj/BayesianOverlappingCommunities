function [logprob,NZQ,suff_stats] = logjoint_opt_2(A,Z,Q,pms,suff_stats,mode,idx,NZQ,M_prev)
% Given a matrix A with its corresponding observation X concept matrix
% Z and feature X concept matrix Q, this function calculates
% the probability of this combination given the parameters.

if nargin < 5
    mode = 0;
    idx = 0;
    NZQ = zeros(size(A,1),size(A,2));
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
    
    % Handle missing data (NaNs)
    A(isnan(A)) = 0;
    suff_stats.NPap(isnan(suff_stats.NPap)) = 0;
    suff_stats.NPam(isnan(suff_stats.NPam)) = 0;
    suff_stats.NPbp(isnan(suff_stats.NPbp)) = 0;
    suff_stats.NPbm(isnan(suff_stats.NPbm)) = 0;
    
    suff_stats.NZp = sum(Z)'+ pms(2,1);
    suff_stats.NZm = no_obs-sum(Z)'+ pms(2,2);
    
    suff_stats.NQp = sum(Q)'+pms(1,1);
    suff_stats.NQm = no_feats-sum(Q)'+pms(1,2);
    
    % Adding the parameters to get
    % final input for beta functions
    
    suff_stats.NPaps = sum(suff_stats.NPap)'+pms(3,1);
    suff_stats.NPams = sum(suff_stats.NPam)'+pms(3,2);
    suff_stats.NPbps = sum(suff_stats.NPbp)'+pms(4,1);
    suff_stats.NPbms = sum(suff_stats.NPbm)'+pms(4,2);
    
elseif mode(1) == 'z'
    %change = str2num(mode(2));
    A(isnan(A)) = 0;
    i = idx(1);
    d = idx(2);
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Z(i,d) is changed.
    %      NZQ(i,:)    +  Z(i,d)*Q(:,d)' - Z_prev(i,d)*Q(:,d)'
    % (previous value) + (change compared to previous value)
    
    Z_diff = Z(i,d)-M_prev(i,d);
    if Z_diff ~= 0
        NZQ(i,:) = NZQ(i,:)+Z_diff*Q(:,d)';
    end
    
    NPap_old = suff_stats.NPap(i,:);
    NPam_old = suff_stats.NPam(i,:);
    NPbp_old = suff_stats.NPbp(i,:);
    NPbm_old = suff_stats.NPbm(i,:);
    
    PZQ =(0<(NZQ(i,:)));
    AZQ = A(i,:).*PZQ;
    
    % #Elements that share concepts and are linked
    suff_stats.NPap(i,:) = AZQ;
    % #Elements that share concepts and aren't linked
    suff_stats.NPam(i,:) = PZQ-AZQ;
    % #Elements that does not share concepts and are linked
    suff_stats.NPbp(i,:) = A(i,:)-AZQ;
    % #Elements that do not share concepts and aren't linked
    suff_stats.NPbm(i,:) = AZQ+1-PZQ-A(i,:); %(1-A(i,:)).*(1-(0<(NZQ(i,:))));
    
    % Find difference from previous state
    suff_stats.NZp(d) = suff_stats.NZp(d)+Z_diff;
    suff_stats.NZm(d) = suff_stats.NZm(d)-Z_diff;
    
    suff_stats.NPaps = suff_stats.NPaps+suff_stats.NPap(i,:)'-NPap_old';
    suff_stats.NPams = suff_stats.NPams+suff_stats.NPam(i,:)'-NPam_old';
    suff_stats.NPbps = suff_stats.NPbps+suff_stats.NPbp(i,:)'-NPbp_old';
    suff_stats.NPbms = suff_stats.NPbms+suff_stats.NPbm(i,:)'-NPbm_old';
elseif mode(1) == 'q'
    % index for changed value of Q
    A(isnan(A)) = 0;
    j = idx(1);
    d = idx(2);
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Q(j,d) is changed.
    %      NZQ(:,j)    +  Z(:,d)*Q(j,d)' - Z(:,d)*Q_prev(j,d)'
    % (previous value) + (change compared to previous value)

    Q_diff=Q(j,d)-M_prev(j,d);
    if Q_diff ~= 0
        NZQ(:,j) = NZQ(:,j)+Z(:,d)*Q_diff;
    end
    
    NPap_old = suff_stats.NPap(:,j);
    NPam_old = suff_stats.NPam(:,j);
    NPbp_old = suff_stats.NPbp(:,j);
    NPbm_old = suff_stats.NPbm(:,j);

    PZQ =(0<(NZQ(:,j)));
    AZQ = A(:,j).*PZQ;
    
    % #Elements that share concepts and are linked
    suff_stats.NPap(:,j) = AZQ;
    % #Elements that share concepts and aren't linked
    suff_stats.NPam(:,j) = PZQ-AZQ;
    % #Elements that does not share concepts and are linked
    suff_stats.NPbp(:,j) = A(:,j)-AZQ;
    % #Elements that do not share concepts and aren't linked
    suff_stats.NPbm(:,j) = AZQ+1-PZQ-A(:,j); %(1-A(i,:)).*(1-(0<(NZQ(i,:))));
    
    
    % Find difference from previous state
    suff_stats.NQp(d) = suff_stats.NQp(d)+Q_diff;
    suff_stats.NQm(d) = suff_stats.NQm(d)-Q_diff;
    
    
    suff_stats.NPaps(j) = suff_stats.NPaps(j)+sum(suff_stats.NPap(:,j)-NPap_old);
    suff_stats.NPams(j) = suff_stats.NPams(j)+sum(suff_stats.NPam(:,j)-NPam_old);
    suff_stats.NPbps(j) = suff_stats.NPbps(j)+sum(suff_stats.NPbp(:,j)-NPbp_old);
    suff_stats.NPbms(j) = suff_stats.NPbms(j)+sum(suff_stats.NPbm(:,j)-NPbm_old);

end

% Calculating the probability

%a1 = betaln(suff_stats.NPaps,suff_stats.NPams)+betaln(suff_stats.NPbps,suff_stats.NPbms);

%a1 = numeric_integration_beta_incomplete(suff_stats.NPaps,suff_stats.NPams,suff_stats.NPbps,suff_stats.NPbms);
a1 = betaintegral_noverbose(suff_stats.NPaps,suff_stats.NPams,suff_stats.NPbps,suff_stats.NPbms,50000);

b1 = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));

a2 = betaln(suff_stats.NZp,suff_stats.NZm)+betaln(suff_stats.NQp,suff_stats.NQm);
b2 = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));



logprob=sum(a1)-no_feats*b1+sum(a2)-no_concepts*b2;

end