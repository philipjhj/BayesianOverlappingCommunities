clear
load('Data/UNESCO/unesco_multi.mat')

collected_unesco = [czeck; japan; korea];
collected_unesco_objs = [N2; N3; N1];
%%
imagesc(collected_unesco);

% % If any columns with only zeros or full
% sum_cls = sum(collected_unesco);
% no_objs = size(collected_unesco,2);
% collected_unesco = collected_unesco(:,0<sum_cls<no_objs);

%%
clear
load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Results\unesco_samples\data_orginal\unesco_full_data.mat')
bs = [0.25 0.25];
rs = [0.25 0.25];

pas = [1000 1]; % pa_plus, pa_minus
pbs = [1 1000]; % pb_plus, pb_minus

pms = [bs; rs; pas; pbs]; %Parameters

tic
[dZ, dQ, logPs] = Gibbs_sampler_opt_2(collected_unesco,10,pms,25);
toc
%%
save('unesco_data_first_analysis','pms','collected_unesco','collected_unesco_objs','dQ','dZ','logPs','Nfeat')

%%
writeFCAFile(collected_unesco,collected_unesco_objs,Nfeat,'unesco_FCA_format.csv')



