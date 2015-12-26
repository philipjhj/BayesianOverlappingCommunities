function [A,Z,Q] = GenerateDataFromModel(obs,feas,concs,bs,rs,pas,pbs)
%GenerateGraph Summary of this function goes here
%   Detailed explanation goes here

% Prio for observations to belong to concepts
Bds = betarnd(bs(1),bs(2),1,concs); 

% Prio for features to belong to concepts
Rds = betarnd(rs(1),rs(2),1,concs); 

% Prio for if an observation and feature sharing a concept has an link 
pa = betarnd(pas(1),pas(2),1,feas);
% Prio for if an observation and feature not sharing a concept has an link
pb = betarnd(pbs(1),pbs(2),1,feas); 

Z = binornd(1,repmat(Bds,obs,1)); % Matrix: Observations X Concepts
Q = binornd(1,repmat(Rds,feas,1)); % Matrix: Features x Concepts

% Generates Graph
A = zeros(obs,feas);
for i = 1:obs
    for j = 1:feas
        ind = (0<Z(i,:)*Q(j,:)'); % indicator function
        A(i,j) = binornd(1,pa(j)^(ind)*pb(j)^(1-ind));
    end
end


end

