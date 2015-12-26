clear

load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Data\worklife_japan_maleB.mat')
japan_obj = subjects_male;
japan_values = VBj; 

load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Data\worklife_norway_maleB.mat')
norway_obj = subjects_male;
norway_values = VBn;

load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Data\worklife_sweden_maleB.mat')
sweden_obj = subjects_male;
sweden_values = VBs;



wvs_japan_norway_sweden_males = [japan_values;norway_values;sweden_values];
obj_names = [japan_obj;norway_obj;sweden_obj];

fea_names = values;

%%

%save('wvs_collected_japan_norway_sweden_males','wvs_japan_norway_sweden_males','fea_names','obj_names')

%%
% 
% %%
% imagesc(wvs_japan_norway_sweden_males);
% 
% % If any columns with only zeros or full
% sum_cls = sum(wvs_japan_norway_sweden_males);
% no_objs = size(wvs_japan_norway_sweden_males,2);
% wvs_japan_norway_sweden_males = wvs_japan_norway_sweden_males(:,0<sum_cls<no_objs);

%%
clear
load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Results\WorldValueSurvey\data_original\wvs_collected_japan_norway_sweden_males.mat')

bs = [0.25 0.25];
rs = [0.25 0.25];

pas = [1000 1]; % pa_plus, pa_minus
pbs = [1 1000]; % pb_plus, pb_minus

pms = [bs; rs; pas; pbs]; %Parameters
tic;
[dZ, dQ, logPs] = Gibbs_sampler_opt_2(wvs_japan_norway_sweden_males,10,pms,25);
toc

%%
save('unesco_data_first_analysis','pms','collected_unesco','collected_unesco_objs','dQ','dZ','logPs','Nfeat')

%%
writeFCAFile(wvs_japan_norway_sweden_males,obj_names,fea_names,'wvs_FCA_format.csv')



