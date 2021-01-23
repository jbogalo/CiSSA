function xe = extend(x,H)
% EXTEND - Extends time series to perform Singular Spectrum Analysis. 
%
% This MATLAB function extends the time series at the beginning and end. 
%
% Syntax:       xe = extend(x,H)
%
% Input arguments:
% x:    Column vector with the original time series.
% H:    A number which determines the extension type.
%
% Output arguments:
% xe:   The extended time series.

% -------------------------------------------------------
% Dimensions
% -------------------------------------------------------
T = length(x);

% -------------------------------------------------------
% Extend
% -------------------------------------------------------
switch H
    case 0      % No extension
        xe = x;
    case T      % Mirroring
        xe = [flipud(x(1:H)); x; flipud(x(end-H+1:end))];
    otherwise   % Autoregressive extension
        % AR coefficients of the differentiated series
        p = fix(T/3);
        dx = diff(x);
        A = aryule(dx,p);
        % Right extension
        y = x;
        dy = diff(y);
        er = filter(A,1,dy);
        dy = filter(1,A,[er; zeros(H,1)]);
        y = y(1)+[0; cumsum(dy)];
        % Left extension
        y = flipud(y);
        dy = diff(y);
        er = filter(A,1,dy);
        dy = filter(1,A,[er; zeros(H,1)]);
        y = y(1)+[0; cumsum(dy)];
        % Extended series
        xe = flipud(y);
end
