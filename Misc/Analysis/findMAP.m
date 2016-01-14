function [max_prob, sample_max_Z, sample_max_Q,i,idx2] = findMAP(logPs)
% Find the sample with highest probability and its value

[value1 idx1]=max(logPs);

[value2 idx2]=max(value1);

max_prob = value2;
i = idx1(idx2);

if idx1(idx2)==1
    sample_max_Q = idx2;
    sample_max_Z = sample_max_Q-1;
else
    sample_max_Z = idx2;
    sample_max_Q = sample_max_Z;
end