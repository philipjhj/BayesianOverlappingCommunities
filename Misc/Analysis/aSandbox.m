%% Intialize
clear
obs = 10; % No. of observations
feas = 10; % No. of features
concs = 5; % No. of concepts

bs = [100 400];
rs = [100 400];

pas = [1000000 0.01]; % pa_plus, pa_minus
pbs = [0.01 1000000]; % pb_plus, pb_minus

[A,Z,Q] = GenerateDataFromModel(obs,feas,concs,bs,rs,pas,pbs);


pms = [bs; rs; pas; pbs;]; %Parameters

%% plot data
figure(2)
[Z_idx,Q_idx]=plotZQ(Z,Q,0,1);
figure(3)
imagesc(A(Z_idx,Q_idx))

%% Create missing data

Percent = 0.05
A(randperm(numel(A),round(numel(A)*Percent))) = NaN




%% Plots
figure(1)
[Z_idx,Q_idx]=plotZQ(Z,Q,0,1);
Z_idx'
Q_idx'

figure(2)
subplot(121)
imagesc(A(Z_idx,Q_idx))
axis image
subplot(122)
imagesc(A)
axis image
%plotAB(A(Z_idx,Q_idx),Z,Q)

%% Joint fordeling


pms = [bs; rs; pas; pbs;]; %Parameters
%%
prob = logjoint(A,Z,Q,pms)

%% Gibbs

pms = [bs; rs; pas; pbs;]; %Parameters

%%

pms = [0.1 0.1; 0.1 0.1; 1 0.1; 0.1 1];
%%

load('davis')
% Initialize hyperparameters
bs = [1 1];
rs = [1 1];
pas = [1 1]; % pa_plus, pa_minus
pbs = [1 1]; % pb_plus, pb_minus
pms = [bs; rs; pas; pbs]; %hyperparameters

T = 1000;
iters = 4;
col = hsv(iters);
%figure(5)
%clf
hold on
test=1;
for i = 1:iters
    tic
    if test == 1
        %[dZ,dQ, logPs{i},inferred_pms{i}] = Gibbs_sampler_opt_2_hyper(A,T,pms,concs);
        %[dZ,dQ, logPs{i}] = Gibbs_sampler_opt_2_2(A,T,pms,4);
        %[dZ,dQ, logPs{i}] = Gibbs_sampler_opt_2_multiflip(A,T,pms,concs);
       [dZ,dQ, logPs{i},inferred_pms{i}] = Gibbs_sampler_final(A,T,pms,2,2);
       %[dZ,dQ, logPs{i},inferred_pms{i}] = Gibbs_cluster(A,T,pms,concs,4);
        %plot(logPs{i}(:,3),'color',col(i,:))
    else
        [dZ,dQ, logPs] = Gibbs_sampler_opt_2(A,T,pms,4);
    end
    i
    toc
end
hold off

%%
figure(5)
for i = 1:4
    hold on
    plot(logPs{i}(:,3),'color',col(i,:))
end

%%
col = 'ymcr';
figure(3)
hold on
maxP = -Inf;
for i = 1:4
    %load(strcat('/media/philipjhj/Data/OneDrive/Studie/0Spec - Bayesian Association Mining/MATLAB/Results_cluster/Results/davis_d6_k3_run',num2str(i),'.mat'))
    load(strcat('/media/philipjhj/Data/OneDrive/Studie/0Spec - Bayesian Association Mining/MATLAB/Results_cluster/Results/AWA_d30_k4_run',num2str(i),'.mat'))
    plot(logPs(:,3),col(i))
    [maxlogP, max_z,max_q]=Find_max(logPs);
    if maxP < maxlogP;
        maxP = maxlogP;
        best_run = i;
        best_z = max_z;
        best_q = max_q;
    end
end

%%

best_run



%%

probs=[];
probs2=[];
for i = 1:T
    probs(i) = logjoint(A,dZ{i},dQ{i},pms)-logPs{end}(i,2);
end
for i = 1:T-1
    probs2(i) = logjoint(A,dZ{i},dQ{i+1},pms)-logPs{end}(i+1,1);
end
figure(4)
clf
plot(probs)
hold on
plot(probs2,'red')
hold off
%% test of conditional vs joint
probs=[]
test=2
figure(5)
if test == 2
        for i = 1:T
            probs(i) = logjoint(A,dZ{i},dQ{i},inferred_pms{i})-logPs(i,3);
        end
else
    for i = 1:T
        probs(i) = logjoint(A,dZ{i},dQ{i},pms)-logPs{1}(i,2);
    end
end

plot(probs)

%%
if test == 2
    figure(3)
else
    figure(4)
end
plot(probs)


%%
for i = 1:1
    [maxlogP(i), max_z(i),max_q(i)]=Find_max(logPs)
end
%%
logjoint(A,dZ{max_z},dQ{max_q},inferred_pms{min([max_q, max_z])})
%%
plot(logPs(:,3))

%%
h=1
figure(1)
plotZQ(Z,Q,0,1)
figure(2)

plotZQ(dZ{max_z(h)},dQ{max_q(h)},0,1)
%plotZQ(dZ{end},dQ{end},1,1)

%%
%h=1
%igure(1)
%plotZQ(Z,Q,1,1)
figure(2)
plotZQ(dZ{max_z},dQ{max_q},0,1)
%plotZQ(dZ{end},dQ{end},1,1)

%%
figure(2)
[Z_idx,Q_idx]=plotZQ(dZ{max_z},dQ{max_q},0,1);
Z_idx';
Q_idx';

figure(1)
subplot(121)
imagesc(A(Z_idx,Q_idx))
axis image
subplot(122)
imagesc(A)
axis image

%%

comm = 1;

women(Z_idx(dZ{max_z}(:,comm)==1))

events(Q_idx(dQ{max_q}(:,comm)==1))



%%
%clear

bs = [0.25 0.25];
rs = [0.25 0.25];

pas = [10 1]; % pa_plus, pa_minus
pbs = [1 10]; % pb_plus, pb_minus

pms = [bs; rs; pas; pbs]; %Parameters
%%
load('C:\Users\Philip\Dropbox\Bachelorprojekt - Philip s113245\MATLAB\Data\worklife_japan_maleB.mat')
%load('ZQ_japan_male_VBj.mat')
A = VBj;
%load('C:\Users\Philip\SkyDrive\Studie\Bachelor Projekt\MATLAB\Results_testing\Gibbs_sampler_testing\Small_data\gibbs_3000samples_logopt_data01_03.mat')
%load('C:\Users\Philip\SkyDrive\Studie\Bachelor Projekt\MATLAB\Results_testing\Gibbs_sampler_testing\Small_data\generated_data01.mat')
%%
[dZ, dQ] = Gibbs_sampler_opt_2(collected_unesco,100,pms,25);

%%
probs = fliplr(rot90(logPs,-1));
%%
probs = reshape(probs,1,size(logPs,1)*size(logPs,2));

%%

plot(probs(1:1:end))


%%
% FIND MAX PROB FROM SAMPLES
max_P = -realmax;
T = length(dZ);

A = wvs_japan_norway_sweden_males;

for i = 1:25:T
    for j = 1:25:T
        P = logjoint(A,dZ{i},dQ{j},pms);
        if  P > max_P
            max_P = P;
            max_i = i;
            max_j = j;
        end
    end
end

%%
dQ_best = dQ{max_j};
dZ_best = dZ{max_i};

max_P
%exp(logjoint(,Z,Q,pms))
[max_i max_j]

%%
subplot(2,2,1)
imagesc(Q)
subplot(2,2,2)
imagesc(dQ_best)

subplot(2,2,3)
imagesc(Z)
subplot(2,2,4)
imagesc(dZ_best)

%% Test k-flip


for k = 2:20%size(Z,2)
    tic
    [a,b,c,d] = k_flip2(Q,Z(1,:),k);
    toc
    if length(b) == 2^k
        fprintf('Correct number of permutations - k = %d\n',k)
    end
    if all(sum(b)==(2^k/2))
        fprintf('Correct permutations - k = %d\n',k)
    end
end