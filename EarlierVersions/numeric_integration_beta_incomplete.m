function [logP2fast]=numeric_integration_beta_incomplete(a1,b1,a0,b0,res)

% function to make numeric integration of the integral:
%  P = \int_0^1 \int_0^\eta1     eta1^(a1-1)*(1-eta1)^(b1-1)*eta0^(a0-1)(1-eta0)^(b0-1) d\eta0 d\eta1
%
% Outputs log(P) calculated using two methods. Method 1 exploits that the
% integral  for \eta0 can be solved using the incomplete beta function.
% Method 2 makes brute force numerical integration on a grid.  Both uses
% res to define the resolution of the grid. (Note that method 2fast is by far the fastest.)
%
% Written by Morten M�rup

if nargin<5
    res=5000;
end

% Method 2 fast:
no_integrations = length(a0);
logQ_tmp=nan(no_integrations,res);
x=linspace(1/(2*res),1-1/(2*res),res);
logP_tmp=((a0-1))*log(x)+((b0-1))*log(1-x);
cumsum_P=zeros(no_integrations,1);
max_old=logP_tmp(:,1);

for t=1:length(logP_tmp)
    max_new=logP_tmp(:,t);                
    max_comp = max_new > max_old;
    
    z = (a1-1)*log(x(t))+(b1-1)*log(1-x(t));
    d1=exp(max_old-max_new);
    d2=exp(max_new-max_old);
    for i=1:no_integrations
        if max_comp(i)
            cumsum_P(i)=cumsum_P(i)*d1(i)+0.5;
            logQ_tmp(i,t)=log(cumsum_P(i))+max_new(i)+z(i); 
            cumsum_P(i)=cumsum_P(i)+0.5;
            max_old(i)=max_new(i);
        else
            q=d2(i);
            cumsum_P(i)=cumsum_P(i)+0.5*q;
            logQ_tmp(i,t)=log(cumsum_P(i))+max_old(i)+z(i);       
            cumsum_P(i)=cumsum_P(i)+0.5*q;
        end
    end
end
logP2fast=log(sum(exp(logQ_tmp-max(logQ_tmp,[],2)*ones(1,res)),2))+max(logQ_tmp,[],2)-2*log(res);
