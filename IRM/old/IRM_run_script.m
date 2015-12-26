
%% IRM
opts.maxiter = 200;
[L,cpu_time,Z,eta,sample,West]=IRMNway_sampler_hyper(collected_unesco,zeros(size(collected_unesco)),[5 5],opts)
%% LogP
figure; plot(L)
%% Mold 1 max + 2
figure, imagesc(sample.MAP.Z{1})
figure, imagesc(sample.MAP.Z{2})
%% Etas
figure, imagesc(sample.MAP.eta); colorbar;
ylabel('Concepts of objects');
xlabel('Concepts of features');
%% Hyperparameters
sample
sample.MAP

%%
%% IRM
opts.maxiter = 200;
[L,cpu_time,Z,eta,sample,West]=IRMNway_sampler_hyper(wvs_japan_norway_sweden_males,zeros(size(wvs_japan_norway_sweden_males)),[5 5],opts)
%% LogP
figure; plot(L)
%% Mold 1 max + 2
figure, imagesc(sample.MAP.Z{1})
figure, imagesc(sample.MAP.Z{2})
%% Etas
figure, imagesc(sample.MAP.eta); colorbar;
ylabel('Concepts of objects');
xlabel('Concepts of features');
%% Hyperparameters
sample
sample.MAP


%%

A = [1 1 1 0 0;
    1 1 1 1 0;
    0 1 1 0 1;
    1 0 1 1 1;
    0 0 0 1 1];

opts.maxiter = 200;
[L,cpu_time,Z,eta,sample,West]=IRMNway_sampler_hyper(A,zeros(size(A)),[5 5],opts)
%% LogP
figure; plot(L)
%% Mold 1 max + 2
figure, imagesc(sample.MAP.Z{1})
figure, imagesc(sample.MAP.Z{2})
%% Etas
figure, imagesc(sample.MAP.eta); colorbar;
ylabel('Concepts of objects');
xlabel('Concepts of features');
%% Hyperparameters
sample
sample.MAP





