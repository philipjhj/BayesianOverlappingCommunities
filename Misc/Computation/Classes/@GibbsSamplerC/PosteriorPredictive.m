function obj = PosteriorPredictive(obj)
S = 20;

[~,jIdx]=ind2sub(size(obj.A),obj.MissingData);
jIdx = unique(jIdx);

[nObs,nFeas] = size(obj.A);

burnin = round(.5*obj.T);
for t = (burnin+1):obj.T
    NZQ = 0<(obj.zSamples{t}*obj.qSamples{t}');
    paSum = zeros(1,nFeas);
    pbSum = zeros(1,nFeas);
    for j = jIdx;%1:nFeas;
        ninp = obj.SuffStatAll{t}.NPaps(j)+obj.pmsSamples{t}(3,1);
        ninm = obj.SuffStatAll{t}.NPams(j)+obj.pmsSamples{t}(3,2);
        noutp = obj.SuffStatAll{t}.NPbps(j)+obj.pmsSamples{t}(4,1);
        noutm = obj.SuffStatAll{t}.NPbms(j)+obj.pmsSamples{t}(4,2);
        
        a=(sum(betaintegral_noverbose(ninp+1,ninm,noutp,noutm)));
        b=(sum(betaintegral_noverbose(ninp,ninm,noutp+1,noutm)));
        c=(sum(betaintegral_noverbose(ninp,ninm,noutp,noutm)));
%         
%        paSum(j) = exp(a)/exp(c);
%        pbSum(j) = exp(b)/exp(c);
       paSum(j) = exp(a-c);
       pbSum(j) = exp(b-c);
%         myJ = find(j==(1:nFeas));
%         for s = 1:S
%             pa = 0;
%             pb = 100;
%             while pa < pb
%                 pa = betarnd(obj.SuffStatAll{t}.NPaps(j)+...
%                     obj.pmsSamples{t}(3,1),...
%                     obj.SuffStatAll{t}.NPams(j)+obj.pmsSamples{t}(3,2));
%                 
%                 pb = betarnd(obj.SuffStatAll{t}.NPbps(j)+...
%                     obj.pmsSamples{t}(4,1),...
%                     obj.SuffStatAll{t}.NPbms(j)+obj.pmsSamples{t}(4,2));
%             end
%             paSum(myJ) = paSum(myJ) + pa;
%             pbSum(myJ) = pbSum(myJ) + pb;
%         end
    end
%     paSum = paSum/S;
%     pbSum = pbSum/S;
    paMat=repmat(paSum,nObs,1);
    pbMat=repmat(pbSum,nObs,1);
%     paMat=repmat(paSum',1,nFeas);
%     pbMat=repmat(pbSum',1,nFeas);
    
    obj.PostPredDistMat(obj.MissingData) = obj.PostPredDistMat(obj.MissingData) + ...
        paMat(obj.MissingData).*NZQ(obj.MissingData);
    
    obj.PostPredDistMat(obj.MissingData) = obj.PostPredDistMat(obj.MissingData) + ...
        pbMat(obj.MissingData).*(1-NZQ(obj.MissingData));
end
obj.PostPredDistMat = obj.PostPredDistMat/(obj.T-burnin);
end