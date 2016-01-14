function prob = logjointNaive(A,Z,Q,pms,cmp)
% Given a matrix A with its corresponding observation X concept matrix
% Z and feature X concept matrix Q, this function calculates
% the probability of this combination given the parameters.

no_obs = size(A,1);
no_feats = size(A,2);
no_concepts_Z = size(Z,2);
no_concepts_Q = size(Q,2);

% NPap = ones(no_feats,1);
% NPam = ones(no_feats,1);
% NPbp = ones(no_feats,1);
% NPbm = ones(no_feats,1);
% 
% NZp = ones(no_concepts_Z,1);
% NZm = ones(no_concepts_Z,1);
% 
% NQp = ones(no_concepts_Q,1);
% NQm = ones(no_concepts_Q,1);


for j = 1:no_feats
    NPap(j) = sum(A(:,j).*any(Z*Q(j,:)',2));
end

for j = 1:no_feats
    NPam(j) = sum((1-A(:,j)).*any(Z*Q(j,:)',2));
end

for j = 1:no_feats
    NPbp(j) = sum(A(:,j).*(1-any(Z*Q(j,:)',2)));
end

for j = 1:no_feats
    NPbm(j) = sum((1-A(:,j)).*(1-any(Z*Q(j,:)',2)));
end

NZp = sum(Z)';
NZm = no_obs-sum(Z)';

NQp = sum(Q)';
NQm = no_feats-sum(Q)';

%     display(NPap)
%     display(NPam)
%     display(NPbp)
%     display(NPbm)
%     display(NZp)
%     display(NZm)
%     display(NQp)
%     display(NQm)
NPap = NPap'+pms(3,1)
NPam = NPam'+pms(3,2)

NPbp = NPbp' + pms(4,1);
NPbm = NPbm' + pms(4,2);

NZp = NZp + pms(1,1);
NZm = NZm + pms(1,2);

NQp = NQp + pms(2,1);
NQm = NQm + pms(2,2);
if nargin > 4
    a11=sum(cmp.NPaps-NPap);
    a12=sum(cmp.NPams-NPam);
    b11=sum(cmp.NPbps-NPbp);
    b12=sum(cmp.NPbms-NPbm);
    q1=sum(cmp.NQp-NQp);
    q2=sum(cmp.NQm-NQm);
    z1=sum(cmp.NZp-NZp);
    z2=sum(cmp.NZm-NZm);
    all = [a11,a12,b11,b12,q1,q2,z1,z2];
    if any(all>10e-6)
        fprintf('Problem with sufficient stat: %d',find(all,1))
        all
        prob = -1;
    else
        %display('all sufficient statistics correct');
        prob = 1;
    end
    return
end

%fprintf('Results:');
%[NPap NPam NPbp NPbm ]
%[NZp NZm NQp NQm]
% Adding the parameters to get 
% final input for beta functions


%[NPap NPam NPbp NPbm]

% Calculating the probability

%prob1 = betaln(NPap,NPam)+betaln(NPbp,NPbm);
%prob1 = numeric_integration_beta_incomplete(NPap,NPam,NPbp,NPbm);
prob1 = betaintegral_noverbose(NPap,NPam,NPbp,NPbm,50000);

frac1 = betaln(pms(3,1),pms(3,2))+betaln(pms(4,1),pms(4,2));

prob2 = betaln(NZp,NZm)+betaln(NQp,NQm);
frac2 = betaln(pms(1,1),pms(1,2))+betaln(pms(2,1),pms(2,2));

prob=sum(prob1)-no_feats*frac1+sum(prob2)-no_concepts_Z*frac2;
end