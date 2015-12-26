function [logP1, logP2, logP2fast]=numeric_integration_beta_incomplete(a1,b1,a0,b0,res)

% function to make numeric integration of the integral:
%  P = \int_0^1 \int_0^\eta1     eta1^(a1-1)*(1-eta1)^(b1-1)*eta0^(a0-1)(1-eta0)^(b0-1) d\eta0 d\eta1
%
% Outputs log(P) calculated using two methods. Method 1 exploits that the
% integral  for \eta0 can be solved using the incomplete beta function.
% Method 2 makes brute force numerical integration on a grid.  Both uses
% res to define the resolution of the grid. (Note that method 2fast is by far the fastest.)
%
% Written by Morten Mørup

if nargin<5
    res=10000;
end

% Method 1:
const=betaln(a0,b0);
t=0;
logP_tmp=nan(1,res);
for x=linspace(1/(2*res),1-1/(2*res),res)
    t=t+1;
    logP_tmp(t)=(a1-1)*log(x)+(b1-1)*log(1-x)+log(betainc(x,a0,b0));
end
logP1=log(sum(exp(logP_tmp-max(logP_tmp))))+max(logP_tmp)+const-log(res);

% Method 2:
t=0;
log2=log(2);
logQ_tmp=nan(1,res);
x=linspace(1/(2*res),1-1/(2*res),res);
logP_tmp_t=(a0-1)*log(x)+(b0-1)*log(1-x);    
for t=1:length(logP_tmp)
    logP_tmp=logP_tmp_t(1:t);
    logP_tmp(end)=logP_tmp(end)-log2;
    logQ_tmp(t)=log(sum(exp(logP_tmp-max(logP_tmp))))+max(logP_tmp)+(a1-1)*log(x(t))+(b1-1)*log(1-x(t));    
end
logP2=log(sum(exp(logQ_tmp-max(logQ_tmp))))+max(logQ_tmp)-2*log(res);

% Method 2 fast:
logQ_tmp=nan(1,res);
x=linspace(1/(2*res),1-1/(2*res),res);
logP_tmp=(a0-1)*log(x)+(b0-1)*log(1-x);    
cumsum_P=0;
max_old=logP_tmp(1);    
for t=1:length(logP_tmp)
    max_new=logP_tmp(t);                
    if max_new>max_old        
        cumsum_P=cumsum_P*exp(max_old-max_new)+0.5;
        logQ_tmp(t)=log(cumsum_P)+max_new+(a1-1)*log(x(t))+(b1-1)*log(1-x(t)); 
        cumsum_P=cumsum_P+0.5;
        max_old=max_new;
    else
        q=exp(max_new-max_old);
        cumsum_P=cumsum_P+0.5*q;
        logQ_tmp(t)=log(cumsum_P)+max_old+(a1-1)*log(x(t))+(b1-1)*log(1-x(t));       
        cumsum_P=cumsum_P+0.5*q;
    end       
end
logP2fast=log(sum(exp(logQ_tmp-max(logQ_tmp))))+max(logQ_tmp)-2*log(res);
