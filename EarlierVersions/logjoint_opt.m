function [logprob,NZQ,suff_stats] = logjoint_opt(A,Z,Q,pms,suff_stats,mode,idx,NZQ,M_prev)
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
   
    suff_stats.NZp = sum(Z)'+ pms(2,1);
    suff_stats.NZm = no_obs-sum(Z)'+ pms(2,2);
    
    suff_stats.NQp = sum(Q)'+pms(1,1);
    suff_stats.NQm = no_feats-sum(Q)'+pms(1,2);

elseif mode(1) == 'z'
    %change = str2num(mode(2));
    i = idx(1);
    d = idx(2);    
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Z(i,d) is changed.
    %      NZQ(i,:)    +  Z(i,d)*Q(:,d)' - Z_prev(i,d)*Q(:,d)' 
    % (previous value) + (change compared to previous value) 
    NZQ(i,:) = NZQ(i,:)+Z(i,d)*Q(:,d)'-M_prev(i,d)*Q(:,d)';
    % todo : udregn Z_diff og check hvis den er forskellig fra nul, ellers
    % skal NZQ ikke ændres : todo
    
    % todo : brug sammenhængen mellem disse udregninger
    
    % #Elements that share concepts and are linked
    suff_stats.NPap(i,:) = A(i,:).*(0<(NZQ(i,:)));
    % #Elements that share concepts and aren't linked
    suff_stats.NPam(i,:) = (1-A(i,:)).*(0<(NZQ(i,:)));
    % #Elements that does not share concepts and are linked
    suff_stats.NPbp(i,:) = A(i,:).*(1-(0<(NZQ(i,:))));
    % #Elements that do not share concepts and aren't linked
    suff_stats.NPbm(i,:) = (1-A(i,:)).*(1-(0<(NZQ(i,:))));
    
    % todo : udregn ændring og tilføj istedet for sum
    suff_stats.NZp = sum(Z)'+pms(1,1);
    suff_stats.NZm = no_obs-sum(Z)'+pms(1,2);
    
elseif mode(1) == 'q'
    % index for changed value of Q
    j = idx(1);
    d = idx(2);
    
    % Element NZQ(i,j) gives the number of shared concepts
    % for observation i and feature j
    % Below is calculated the change of the number of shared concepts
    % when Q(j,d) is changed.
    %      NZQ(:,j)    +  Z(:,d)*Q(j,d)' - Z(:,d)*Q_prev(j,d)' 
    % (previous value) + (change compared to previous value)
    
    % todo : udregn Z_diff og check hvis den er forskellig fra nul, ellers
    % skal NZQ ikke ændres : todo
    NZQ(:,j) = NZQ(:,j)+Z(:,d)*Q(j,d)'-Z(:,d)*M_prev(j,d)';
    
    
    
    % todo : brug sammenhængen mellem disse udregninger
    % #Elements that share concepts and are linked
    suff_stats.NPap(:,j) = A(:,j).*(0<(NZQ(:,j)));
    % #Elements that share concepts and aren't linked
    suff_stats.NPam(:,j) = (1-A(:,j)).*(0<(NZQ(:,j)));
    % #Elements that does not share concepts and are linked
    suff_stats.NPbp(:,j) = A(:,j).*(1-(0<(NZQ(:,j))));
    % #Elements that do not share concepts and aren't linked
    suff_stats.NPbm(:,j) = (1-A(:,j)).*(1-(0<(NZQ(:,j))));
    
    % todo : udregn ændring og tilføj istedet for sum
    suff_stats.NQp = sum(Q)'+pms(2,1);
    suff_stats.NQm = no_feats-sum(Q)'+pms(2,2);

end

% Adding the parameters to get
% final input for beta functions
NPaps = sum(suff_stats.NPap)'+pms(3,1);
NPams = sum(suff_stats.NPam)'+pms(3,2);
% todo : udregn ændring fra forrige udregning
NPbps = sum(suff_stats.NPbp)' + pms(4,1);
NPbms = sum(suff_stats.NPbm)' + pms(4,2);




% Calculating the probability

a1 = betaln(NPaps,NPams)+betaln(NPbps,NPbms);
b1 = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));

a2 = betaln(suff_stats.NZp,suff_stats.NZm)+betaln(suff_stats.NQp,suff_stats.NQm);
b2 = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));

logprob=sum(a1)-no_feats*b1+sum(a2)-no_concepts*b2;

end