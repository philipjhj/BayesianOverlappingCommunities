% Sript to load data and time the gibbs sampler
clear all;
load Results_testing\Gibbs_sampler_testing\generated_200_200_100_data01.mat
d = size(Z,2);
T = 1;

profile on
[dZ,dQ] = Gibbs_sampler_opt(A,T,pms,d);
profile off

c = clock;
time = strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)));

filename = strcat('gibbs_',num2str(T),'_samples_logopt_data_200_200_100_',time)
profilename = strcat('p_',num2str(T),'_samples_logopt_data_200_200_100_',time)

p = profile('info');

save(profilename,'p');
save(filename,'dZ','dQ');