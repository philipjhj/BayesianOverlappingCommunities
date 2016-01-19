% IRM script


% LEAVE-ONE-OUT
%load('C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianassociationmining\Data\davis\davis.mat')
A = gibbs.A;
PostPred = zeros(size(A));

for i = 1:size(A,1);
    for j = 1:size(A,2);
        disp([i j])
        W = zeros(size(A));
        W(i,j) = 1;
        [L,~,Z,eta,sample,West]=IRMNway(A,W,[2 2]);
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

% leave 0.x out

W = zeros(size(testclassL.A));
W(testclassL.MissingData) = 1;

[L,~,Z,eta,sample,West]=IRMNway(testclassL.Atrue,W,[5 5]);
PostPred = West;
md = testclassL.MissingData;
[X,Y,T,AUC] = perfcurve(reshape(testclass.A(md),1,numel(testclass.A(md))),reshape(PostPred(md),1,numel(PostPred(md))),1);
AUC