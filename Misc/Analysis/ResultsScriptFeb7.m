
datafiles = {'/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestDataFromModelPerformanceNoisy_d20_k3_run',...
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestDataFromModelPerformanceStrong_d20_k3_run',...
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestNumberOfComponentsNoisy_d50_k3_run',...
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestNumberOfComponentsStrong_d50_k3_run',...
    'AWA_d12_k3_run',...
    'davis_d4_k3_run',...
    'MouseBrain_contra_d12_k3_run',...
    'MouseBrain_ipsi_d12_k3_run'};

all_AUC = zeros(10,3);
%%
for i = 1:10
    disp(i)
    load(strcat(datafiles{5},int2str(i),'.mat'))
    %load(strcat('/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/davis_d5_k3_run',int2str(i),'.mat'))
    
    idx=mod(i-1,10)+1;
    
    gibbs=gibbs.PosteriorPredictive;
    gibbs=gibbs.computeAUC;
    all_AUC(idx,ceil(i/10)) = gibbs.AUC;
    
    
    for j = 1:5
        [irm,bern]=calc_AUC_IRMBERNMIX(gibbs);
            
        all_AUC(idx,2)=all_AUC(idx,2)+irm;
        all_AUC(idx,3)=all_AUC(idx,3)+bern;
    end
    all_AUC(idx,2)=all_AUC(idx,2)/5;
    all_AUC(idx,3)=all_AUC(idx,3)/5;
end

%%
clf;
means= mean(all_AUC);
means= repmat(means,2,1);
figure(4)
b=bar(all_AUC);

b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'y';
hold on
plot(repmat([0.5; 10.5],1,3),means,'LineWidth',3)

ylim([(min(min(all_AUC))-0.1) (max(max(all_AUC))+0.1)])
legend('BAM','IRM','BERNMIX')
xlabel('Data Sets')
ylabel('AUC')


%%

for i = 1:10
    disp(i)
    load(strcat('/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestNumberOfComponentsNoisy_d50_k3_run',int2str(i),'.mat'))
    gibbs.plotZQ
    pause(2)
end

%%


clear


datafiles = {'/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestDataFromModelPerformanceNoisy_d20_k3_run',
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestDataFromModelPerformanceStrong_d20_k3_run',
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestNumberOfComponentsNoisy_d50_k3_run',
    '/media/data/Dropbox/Studie/Bayesian Associaion Model - Philip/code/Output/TestNumberOfComponentsStrong_d50_k3_run'};
for j = 3
    for i = 21:29
        load(strcat(datafiles{j},int2str(i),'.mat'))
        
        disp(i)
        assert(exist('gibbs','var')==1,'Gibbs does not exist.')
        
        gibbs = gibbs.computeZQMAP;
        figure(3)
        gibbs.plotASorted
        figure(1)
        gibbs.plotZQ
        figure(2)
        gibbs.plotChains
        pause(2)
    end
end
%%




