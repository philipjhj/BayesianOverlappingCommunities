% Script to generate a specific data size and time gibbs sampler
clear all;
obs = 500; % No. of observations
feas = 500; % No. of features
concs = 25; % No. of concepts

bs = [0.25 0.25];
rs = [0.25 0.25];

pas = [10 1]; % pa_plus, pa_minus
pbs = [1 10]; % pb_plus, pb_minus

[A,Z,Q] = GenerateGraphSimple(obs,feas,concs,bs,rs,pas,pbs);

pms = [bs; rs; pas; pbs;]; %Parameters

d = concs;
T = 1;
profile on
[dZ,dQ] = Gibbs_sampler_opt(A,T,pms,d);
profile off


%c = clock;
%time = strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)));
filename = strcat('gibbs_',num2str(T),'_concepts_logopt_data_',num2str(obs),'_',num2str(feas),'_',num2str(concs))
profilename = strcat('p_',num2str(T),'_concepts_logopt_data_',num2str(obs),'_',num2str(feas),'_',num2str(concs))c

p = profile('info');

save(profilename,'p');

save(filename,'dZ','dQ');