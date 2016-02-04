% This script contains code to generated results

%% Davis *********
for i = 3
  %  load(strcat(['C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianas'...
   %     'sociationmining\Output\davis_d5_k3_run'],int2str(i),'.mat'));
     load(strcat(['C:\Users\Philip H. Jorgensen\Documents\MATLAB\bayensianas'...
        'sociationmining\Output\AWA_d30_k3_run'],int2str(i),'.mat'));
%     load(strcat(['C:\Users\Philip H. Jorgensen\Documents\MATLAB\'...
%     'bayensianassociationmining\Output\drugs_d20_k3_run'],int2str(i),'.mat'));
% 
%     gibbs.zLabels = women;
%     gibbs.qLabels = events;

    gibbs.zLabels = animals;
    gibbs.qLabels = predicates;
    
%     gibbs.zLabels = drug_list;
%     gibbs.qLabels = se_list;

    gibbs = gibbs.PosteriorPredictive;
    gibbs = gibbs.computeAUC;
    gibbs.AUC

    figure(1)
    gibbs.plotASorted
    figure(2)
    gibbs.CompareTrueSampleZQ(1,2)
    figure(3)
    gibbs.CompareTrueSampleZQ(0,2)
    figure(4)
    clf
    try
        gibbs.plotChains
    catch
        
    end
    pause(2)
end


%% Synthetic **********
% Dette er godt eksemple paa at modellen virker
% config 20x20x3, 2000 samples, rng 6
% bs = [10 10];
% rs = [10 10];
% pas = [1 0.1]; % pa_plus, pa_minus
% pbs = [0.1 1]; % pb_plus, pb_minus

myrng=6;
gibbs = GibbsSamplerC(myrng)
figure(4)
subplot(2,2,[1 3])
imagesc(gibbs.Atrue)
subplot(2,2,2)
imagesc(gibbs.zTrue)
subplot(2,2,4)
imagesc(gibbs.qTrue)
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
    gibbs.CompareTrueSampleZQ(1,2)
    figure(3)
    clf
    gibbs.plotChains
end