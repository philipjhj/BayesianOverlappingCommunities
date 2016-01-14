%% Intialize
obs = 2; % No. of observations
feas = 2; % No. of features
concs = 2; % No. of concepts

bs = [0.35 0.35];
rs = [0.35 0.35];

pas = [10 1]; % pa_plus, pa_minus
pbs = [1 10]; % pb_plus, pb_minus


pms = [bs; rs; pas; pbs;]; %Parameters


% the data (slightly different)
v={
    0:1
    0:1
    0:1
    0:1
    };
% the engine
n=numel(v);
x=cell(n,1);
[x{1:n,1}]=ndgrid(v{end:-1:1});
p=reshape(cat(n+1,x{:}),[],n);
p=p(:,end:-1:1);
% the result
% disp(p);

max_iter = size(p,1);
n = sqrt(size(p,2));
%prob = 0;

for i = 1:max_iter
    A = reshape(p(i,:),n,n);
    for j = 1:max_iter
        Z = reshape(p(j,:),n,n);
        for k = 1:max_iter
            Q = reshape(p(k,:),n,n);
            prob(i,j,k) = exp(logjointNaive(A,Z,Q,pms));
        end
    end
end
%prob
jointSum = sum(sum(sum(prob)));
if 10e-03 > 1-jointSum
    display('The sum of the joint is one and therefore correct')
else
    display(jointSum)
end


% %%
% 
% prob_vec = reshape(prob,1,max_iter^3)
% 
% plot(cumsum(prob_vec))
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% figure
% mySpyPlot(A);
% 
% %%
% 
% 
