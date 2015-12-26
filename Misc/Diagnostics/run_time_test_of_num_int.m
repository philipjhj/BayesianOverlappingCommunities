T = 1;

%for i=1:10
tic
[dZ,dQ, logPs] = Gibbs_sampler_opt_2_2(A,T,pms,4);
toc
%logPs
%end

%%
% 
% for i=1:80
tic
numeric_integration_beta_incomplete(ones(100,1)*2,ones(100,1)*2,ones(100,1)*2,ones(100,1)*2)
toc
% end