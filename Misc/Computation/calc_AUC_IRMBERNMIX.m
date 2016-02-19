function [IRM_AUC, BERNMIX_AUC] = calc_AUC_IRMBERNMIX(gibbs)


W = zeros(size(gibbs.A));
W(gibbs.MissingData) = 1;

[L,~,Z,eta,sample,West]=IRMNway(gibbs.Atrue,W,[20 20]);
PostPred = West;
md = gibbs.MissingData;

[X,Y,T,IRM_AUC] = perfcurve(reshape(gibbs.Atrue(md),1,numel(gibbs.Atrue(md))),reshape(PostPred(md),1,numel(PostPred(md))),1);

[L,~,Z,eta,sample,West]=IRMNway_Philip(gibbs.Atrue,W,[20 20]);
PostPred = West;
md = gibbs.MissingData;
[X,Y,T,BERNMIX_AUC] = perfcurve(reshape(gibbs.Atrue(md),1,numel(gibbs.Atrue(md))),reshape(PostPred(md),1,numel(PostPred(md))),1);

end
