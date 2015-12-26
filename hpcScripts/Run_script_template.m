% Input
% k=4;
% d = 50;
% T = 3;
% nFile = 2;
% method = 'kflip';
% save_f = 1;
% iter = job_no (enviroment variable)
% missing_data = 1;
% Percent = 0.10;

% Names for input and output files
datafile = {'drugsandsideeffects_data.mat',
    'davis.mat',
    'AWA_data.mat',
    'synthetic_data.mat'};
filename_prefix = {'drugs','davis','AWA','synthetic'};

% Load Data
addpath(genpath('Functions'),genpath('Data'))
load(datafile{nFile});

if missing_data == 1
    %Percent = 0.05
    A_orig = A;
    A(randperm(numel(A),round(numel(A)*Percent))) = NaN;
else
    
    
end


% Initialize hyperparameters
bs = [1 1];
rs = [1 1];
pas = [1 1]; % pa_plus, pa_minus
pbs = [1 1]; % pb_plus, pb_minus
pms = [bs; rs; pas; pbs]; %hyperparameters

% Run the sampler
tic;
if strcmp(method,'kflip')
    [dZ, dQ, logPs, inferred_pms,suffstats,cpu_time,wall_time] = Gibbs_sampler_final(A,T,pms,d,k);
elseif strcmp(method,'normgibbs')
    [dZ, dQ, logPs, inferred_pms] = Gibbs_sampler_opt_2_hyper(A,T,pms,d);
end
totalruntime=toc;

% Save output
if save_f
    filename = strcat('Results//',filename_prefix{nFile},'_d',num2str(d),'_k',num2str(k),'_run',num2str(iter));
    save(filename);
end




