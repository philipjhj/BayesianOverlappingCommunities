function plotAB(A,Z,Q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
figure(6)
subplot(3,1,2)
imagesc(A)

title('Graph')

% Observations X Features Matrix
B = 0<Z*Q';

subplot(3,1,1)
imagesc(B)
title('Observation X Features')


subplot(3,1,3)
imagesc(abs(A-B))

end

