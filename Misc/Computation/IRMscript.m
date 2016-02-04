% IRM script


% LEAVE-ONE-OUT
load('C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianassociationmining\Data\davis\davis.mat')
%A = gibbs.A;
PostPred = zeros(size(A));

for i = 1:size(A,1);
    for j = 1:size(A,2);
        disp([i j])
        W = zeros(size(A));
        W(i,j) = 1;
        [L,~,Z,eta,sample,West]=IRMNway_Philip(A,W,[2 2]);
        PostPred = PostPred+West;
    end
end

%%
% AUC

[X,Y,T,AUC] = perfcurve(reshape(A,1,numel(A)),reshape(PostPred,1,numel(PostPred)),1);
AUC

%figure(3)
%plot(X,Y)

%%
% Plot matrics to compare

figure(1)
imagesc(PostPred)
colorbar
figure(2)
imagesc(A)

%% leave 0.x out
    
W = zeros(size(gibbs.A));
W(gibbs.MissingData) = 1;

[L,~,Z,eta,sample,West]=IRMNway_Philip(gibbs.Atrue,W,[20 20]);
PostPred = West;
md = gibbs.MissingData;
[X,Y,T,AUC] = perfcurve(reshape(gibbs.Atrue(md),1,numel(gibbs.Atrue(md))),reshape(PostPred(md),1,numel(PostPred(md))),1);
AUC