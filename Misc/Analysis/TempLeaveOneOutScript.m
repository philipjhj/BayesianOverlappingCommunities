load('C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianassociationmining\Data\davis\davis.mat')

PostPred = zeros(size(A));

for i = 1:18;
    for j = 1:14;
        testclass = GibbsSamplerC(A,4,2);
        testclass.verbose = false;
        testclass=testclass.createMissingdata(i,j);
        testclass=testclass.Sampler(100);
        testclass=testclass.PosteriorPredictive;
        PostPred = PostPred+testclass.PostPredDistMat;
        display(j)
    end
    display(i)
end
PostPred
%%
figure(1)
imagesc(PostPred)
colorbar
figure(2)
imagesc(testclass.A)
%%

[X,Y,T,AUC] = perfcurve(reshape(testclass.A,1,numel(testclass.A)),reshape(PostPred,1,numel(PostPred)),1);
AUC
figure(3)
plot(X,Y)