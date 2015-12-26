% BETAINTEGRAL Evaluates the integral
% $\int_0^1 \int_0^y y^{a-1}(1-y)^{b-1}x^{c-1}(1-x)^{d-1}\mathrm{d}y\mathrm{d}x$
%
% logP = betaintegral(a, b, c, d) evaluates the integral.
%
% logP = betaintegral(a, b, c, d, imax) evaluates the integral using no
% more than imax terms. If imax is not specified it is set to 10000
%
% Implemented in betaintegral.cpp. Compile for MATLAB using the command
%   mex betaintegral.cpp
%
% Copyright 2015, Mikkel N. Schmidt, mnsc@dtu.dk