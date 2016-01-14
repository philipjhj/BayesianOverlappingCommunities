%% Intialize
clear
obs = 4; % No. of observations
feas = 4; % No. of features
concs = 2; % No. of concepts

bs = [1 1];
rs = [1 1];

pas = [1 0.1]; % pa_plus, pa_minus
pbs = [0.1 1]; % pb_plus, pb_minus
rng(2)
[A,zTrue,qTrue] = GenerateDataFromModel(obs,feas,concs,bs,rs,pas,pbs);


pms = [bs; rs; pas; pbs;];%Parameters
%
%pmsEstimate = [50 100; 50 100; 1000 0.1; 0.1 1000];
%
T = 400;
myGibbsChain = GibbsChainData(A,T,concs,2,pms,zTrue,qTrue);

nA = numel(A);
myResults = cell(1,nA);


activeGibbs=true;
TestError = true;
ConstrainedInt = true;
activeHyper = true;

reconstructA = zeros(size(A));

for i = 1:nA
    
    tempGibbsChain = myGibbsChain;
    tempGibbsChain.A(i) = nan;
    
    myResults{i} = GibbsSampler(tempGibbsChain,activeGibbs,TestError,ConstrainedInt,activeHyper);
    display(i)
    logp0 = 0;
    logp1 = 0;
    burnin=100;
    usedT = T-burnin;
    P = zeros(1,usedT);
    for j = (burnin+1):T
        suff_stats_c = SuffStatsClass;
        Azero = myGibbsChain.A;
        Azero(i) = 0;
        logp0 = LogJoint(ConstrainedInt,Azero,myResults{i}.zSamples{j},myResults{i}.qSamples{j},myResults{i}.pms,suff_stats_c);
        logp1 = LogJoint(ConstrainedInt,myGibbsChain.A,myResults{i}.zSamples{j},myResults{i}.qSamples{j},myResults{i}.pms,suff_stats_c);
        
        ps=[logp0 logp1];
        Ps = exp(ps-max(ps));
        P(j) = 1-Ps(1)/sum(Ps);
    end
    
    reconstructA(i) = mean(P); 
end

%%
figure(3)
imagesc(A)
figure(1)
imagesc(reconstructA>0.5,[0 1])
colorbar

%%

%plot(myResults{1}.logPs(1,:))
for i = 1:nA
    figure(2)
    imagesc(myResults{i}.A)
    title(sprintf('plot %d',i))
    pause(0.5)
end
%%
figure(3)
imagesc(myGibbsChain.A)
