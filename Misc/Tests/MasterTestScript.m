% TEST SAMPLER on leave-one-out
% Config: (10x10x3, rng(1), 100 samples, 0.5 burnin)
% gibbs = GibbsSamplerC(2)
load('C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianassociationmining\Data\davis\davis.mat')

gibbs = GibbsSamplerC(A,5,3)
gibbs.verbose = 0
gibbs=gibbs.LeaveOneOut(50)

gibbs=gibbs.computeAUC;
gibbs.AUC
%% Test generation of data
% for i = 1:20
myrng=2;
gibbs = GibbsSamplerC(myrng)
figure(4)
subplot(2,2,[1 3])
imagesc(gibbs.Atrue)
title(sprintf('rng: %d',myrng))
subplot(2,2,2)
imagesc(gibbs.zTrue)
subplot(2,2,4)
imagesc(gibbs.qTrue)
pause(1)
end
%% Test long run
continueGibbs = 0;

nSamples = 2000;

for i = 1:1
    disp(i)
    if ~continueGibbs
        gibbs = GibbsSamplerC(myrng)
        gibbs.Debug = 1;
        gibbs.verbose = 1;
    end
    gibbs = gibbs.Sampler(nSamples)
    figure(1)
    clf
    gibbs.plotASorted
    figure(2)
    clf
    gibbs.CompareTrueSampleZQ(0,2)
    figure(3)
    clf
    gibbs.plotChains
end

%%
figure(1)
clf
gibbs.CompareTrueSampleZQ

