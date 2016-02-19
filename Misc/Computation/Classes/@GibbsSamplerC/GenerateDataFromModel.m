function [A,Z,Q,pa,pb,Bds,Rds] = GenerateDataFromModel(obs,feas,concs,bs,rs,pas,pbs)
%GenerateGraph Summary of this function goes here
%   Detailed explanation goes here

% Prio for observations to belong to concepts
Bds = betarnd(bs(1),bs(2),1,concs);

% Prio for features to belong to concepts
Rds = betarnd(rs(1),rs(2),1,concs);
pa = -1;
pb = 0;

while pa < pb
    % Prio for if an observation and feature sharing a concept has an link
    pa = betarnd(pas(1),pas(2),1,feas);
    % Prio for if an observation and feature not sharing a concept has an link
    pb = betarnd(pbs(1),pbs(2),1,feas);
end
Z = binornd(1,repmat(Bds,obs,1)); % Matrix: Observations X Concepts
Q = binornd(1,repmat(Rds,feas,1)); % Matrix: Features x Concepts

% Generates Graph
%A = zeros(obs,feas);
NZQ = (0<Z*Q');
p = bsxfun(@times,NZQ,pa)+bsxfun(@times,1-NZQ,pa); % indicator function
A = binornd(1,p);


end

