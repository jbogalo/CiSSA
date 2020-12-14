function [Z, psd] = cissa(x,L,varargin)
% CiSSA - Circulant Singular Spectrum Analysis.
%
% This MATLAB function returns the elementary reconstructed series by
% frequency, Z, and the estimated power spectral density, psd, of the
% input time series, x, using Circulant Singular Spectrum Analysis
% given a window length, L.
%
% Syntax:     [Z, psd] = cissa(x,L)
%             [Z, psd] = cissa(x,L,H)
%
% Input arguments:
% x: Columm vector containing the time series original data.
% L: Window length.
% H: Optional. Related to the characteristics of the time series.
%    H=0 Autoregressive extension (default). It is indicated for stationary
%        and stochastic trend time series as well.
%    H=1 Mirroing. It can be used with stationary time series and works
%        well for AM-FM signals.
%    H=2 No extension. It is suitable for deterministic time series.
%
% Output arguments:
% Z:   Matrix whose columns are the reconstructed components by frequency.
% psd: Column vector with the estimated power spectral density at
%      frequencies w(k)=(k-1)/L, k=1,2,...,L. This is, the eigenvalues of
%      the circulant matrix of second moments.
%
% See also: group
%
% -------------------------------------------------------------------------
% References:
% [1] Bógalo, J., Poncela, P., and Senra, E. "Circulant Singular Spectrum
%     Analysis: A new automated procedure for signal extraction". Signal
%     Processing. Vol. 179, 2021, in progress.
%     https://doi.org/10.1016/j.sigpro.2020.107824.
% -------------------------------------------------------------------------

% -------------------------------------------------------
% Checking the input arguments
% -------------------------------------------------------
% Dimensions
if size(x,1)==1
    x = x';
end
T = length(x);
N = T-L+1;
if L>N
    error('***  The window length must be less than T/2  ***');
end

% Type of extension depending on H
if nargin>2
    switch varargin{1}
        case 1
            H = T;
        case 2
            H = 0;
        otherwise
            H = L;
    end
else
    H = L;
end

% Number of symmetryc frequency pairs around 1/2
if mod(L,2)
    nf2 = (L+1)/2-1;
else
    nf2 = L/2-1;
end

% Number of frequencies <= 1/2
nft = nf2+abs(mod(L,2)-2);

% -------------------------------------------------------
% Trajectory matrix
% -------------------------------------------------------
% Extended series
xe = extend(x,H);

% Trajectory matrix
col = xe(1:L);
row = xe(L:end);
X = hankel(col,row);
clear col row

% -------------------------------------------------------
% Decomposition
% -------------------------------------------------------
% Autocovariance function
gam = zeros(L,1);
for k=0:L-1
    gam(k+1) = (x(1:T-k)-mean(x))'*(x(1+k:T)-mean(x))/(T-k);
end

% Symmetric Toeplitz covariance matrix S and equivalent circulant matrix C
S = gam(1)*eye(L); C = S;
for i=1:L
    for j=i+1:L
        k = abs(i-j);
        S(i,j) = gam(k+1); S(j,i) = S(i,j);
        C(i,j) = ((L-k)/L)*gam(k+1)+(k/L)*gam(L-k+1); % Pearl (1973)
        C(j,i) = C(i,j);
    end
end
clear gam

% Eigenvectors of circulant matrix (unitary base)
U = dftmtx(L)/sqrt(L);

% Real eigenvectors (orthonormal base)
U(:,1) = real(U(:,1));
for k=1:nf2
    u_k = U(:,k+1);
    U(:,k+1) = sqrt(2)*real(u_k);
    U(:,L+2-(k+1)) = sqrt(2)*imag(u_k);
end
if ~mod(L,2)
    U(:,nft) = real(U(:,nft));
end

% Eigenvalues of circulant matrix: estimated power spectral density
psd = abs(diag(U'*C*U));

% Principal components
W = U'*X;

% -------------------------------------------------------
% Reconstruction
% -------------------------------------------------------
% Elementary reconstructed series
R = zeros(T+2*H,L);
for k=1:L
    R(:,k) = diagaver(U(:,k)*W(k,:));
end

% -------------------------------------------------------
% Grouping by frequency
% -------------------------------------------------------
% Elementary reconstructed series by frequency
Z = zeros(T+2*H,nft);
Z(:,1) = R(:,1);
for k=1:nf2
    Z(:,k+1) = R(:,k+1)+R(:,L+2-(k+1));
end
if ~mod(L,2)
    Z(:,nft) = R(:,nft);
end
Z = Z(H+1:end-H,:);
