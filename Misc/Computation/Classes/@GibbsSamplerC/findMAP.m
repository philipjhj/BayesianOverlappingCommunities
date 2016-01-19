function [max_prob, sample_max_Z, sample_max_Q,i,j] = findMAP(obj)
% Find the sample with highest probability and its value

[value1, idx1]=max(obj.logPs);

[value2, idx2]=max(value1);

max_prob = value2;
j = idx2;
i = idx1(j);

if i==1
    sample_max_Q = j;
    sample_max_Z = sample_max_Q-1;
else
    sample_max_Z = j;
    sample_max_Q = sample_max_Z;
end