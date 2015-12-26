classdef SuffStatsClass
    properties
        NPap
        NPam
        NPbp
        NPbm
        NZp
        NZm
        NQp
        NQm
        % sums of the above ( still pr. feature)
        NPaps
        NPams
        NPbps
        NPbms
        % log-beta evaluated
        betahyper_ZQ
        betahyper_probs
        betasuff_AZQ
        betasuff_sums
    end
end