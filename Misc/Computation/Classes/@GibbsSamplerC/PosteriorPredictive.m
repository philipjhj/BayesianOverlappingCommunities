function obj = PosteriorPredictive(obj)
S = 20;
jIdx = unique(obj.MissingData(:,2));
iIdx = unique(obj.MissingData(:,1));

nFeas = numel(jIdx);
nObs = numel(iIdx);

burnin = 50;
for t = (burnin+1):obj.T
    NZQ = 0<(obj.zSamples{t}*obj.qSamples{t}');
    paSum = zeros(nFeas,1);
    pbSum = zeros(nFeas,1);
    for j = jIdx;
        myJ = find(j==jIdx);
        for s = 1:S
            pa = 0;
            pb = 100;
            while pa < pb
                pa = betarnd(obj.SuffStatAll{t}.NPaps(j)+...
                    obj.pmsSamples{t}(3,1),...
                    obj.SuffStatAll{t}.NPams(j)+obj.pmsSamples{t}(3,2));
                
                pb = betarnd(obj.SuffStatAll{t}.NPbps(j)+...
                    obj.pmsSamples{t}(4,1),...
                    obj.SuffStatAll{t}.NPbms(j)+obj.pmsSamples{t}(4,2));
            end
            paSum(myJ) = paSum(myJ) + pa;
            pbSum(myJ) = pbSum(myJ) + pb;
        end
    end
    paSum = paSum/S;
    pbSum = pbSum/S;
    
    obj.PostPredDistMat(iIdx,jIdx) = obj.PostPredDistMat(iIdx,jIdx) + ...
        repmat(paSum,1,nObs).*NZQ(iIdx,jIdx);
    
    obj.PostPredDistMat(iIdx,jIdx) = obj.PostPredDistMat(iIdx,jIdx) + ...
        repmat(pbSum,1,nObs).*(1-NZQ(iIdx,jIdx));
end
obj.PostPredDistMat = obj.PostPredDistMat/(obj.T-burnin);
end