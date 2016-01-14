function reconstructionA = predictiveDistribution(A,Z,Q,hyperpms,suff_stats)

T=numel(Z);
nNan = numel(find(isnan(A)));
prob_aij = zeros(1,nNan);

for h = find(isnan(A))
    for i = 1:T
        prob_aij(h) =+ exp(sum(betaintegral_noverbose(suff_stats.NPaps+1,suff_stats.NPams,suff_stats.NPbps,suff_stats.NPbms,50000))-no_feats*suff_stats.betahyper_probs);  
    end
    prob_aij(h) = prob_aij(h)/T;
end

A(find(isnan(A))) = prob_aij;
reconstructionA=A;
end