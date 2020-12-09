function xe = extend(x,K)
% EXTEND - Extends time series to perform Singular Spectrum Analysis. 
%
% This MATLAB function extends the time series at the beginning and end. 
%
% Syntax:       xe = extend(x,K)
%
% Input arguments:
% x:    Column vector with the original time series.
% K:    A number which determines the extension type.
%
% Output arguments:
% xe:   The extended time series.

% -------------------------------------------------------
% Dimensions
% -------------------------------------------------------
T = length(x);
p = fix(T/2);

% -------------------------------------------------------
% Extend
% -------------------------------------------------------
switch K
    case 0      % Deterministic
        xe = x;
    case T      % AM-FM signal
        xe = [flipud(x(1:K)); x; flipud(x(end-K+1:end))];
    otherwise   % Time series with a stochastic trend
        % AR coefficients of the differentiated series
        dx = diff(x);
        A = aryule(dx,p);
        % Right extension
        y = x;
        dy = diff(y);
        er = filter(A,1,dy);
        dy = filter(1,A,[er; zeros(K,1)]);
        y = y(1)+[0; cumsum(dy)];
        % Left extension
        y = flipud(y);
        dy = diff(y);
        er = filter(A,1,dy);
        dy = filter(1,A,[er; zeros(K,1)]);
        y = y(1)+[0; cumsum(dy)];
        % Extended series
        xe = flipud(y);
end
