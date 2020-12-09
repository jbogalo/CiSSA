function y = diagaver(Y)
% DIAGAVER - Diagonal averaging for Singular Spectrum Analysis.
%
% This MATLAB function transforms the matrix, Y, into the time series,
% y, by diagonal averaging. This entails averaging the elements of Y over
% its antidiagonals.
%
% Syntax:     y = diagaver(Y)
%
% Input arguments:
% Y:    Matrix associated with a group.
%
% Output arguments:
% y:    Column vector with the reconstructed component.

% -------------------------------------------------------
% Realignment
% -------------------------------------------------------
[LL, NN] = size(Y);
if LL>NN
    Y = Y';
end

% -------------------------------------------------------
% Dimensions
% -------------------------------------------------------
L = min(LL,NN);
N = max(LL,NN);
T = N+L-1;

% -------------------------------------------------------
% Diagonal averaging
% -------------------------------------------------------
y = zeros(T,1);
for t=1:T
    if 1<=t && t<=L-1
        j_inf = 1; j_sup = t;
    elseif L<=t && t<=N
        j_inf = 1; j_sup = L;
    else
        j_inf = t-N+1; j_sup = T-N+1;
    end
    nsum = j_sup-j_inf+1;
    for m=j_inf:j_sup
        y(t) = y(t)+Y(m,t-m+1)/nsum;
    end
end
