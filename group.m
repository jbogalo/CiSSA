function [rc, sh, kg] = group(Z,psd,I)
% GROUP - Grouping step of CiSSA.
%
% This MATLAB function groups the reconstructed components by frequency
% obtained with CiSSA into disjoint subsets and computes the share of the
% corresponding PSD.
%
% Syntax:     [rc, sh, kg] = group(Z,psd,I)
%
% Input arguments:
% Z:   Matrix whose columns are the reconstructed components by frequency
%      obtained with CiSSA.
% psd: Column vector with the estimated power spectral density at
%      frequencies w(k)=(k-1)/L, k=1,2,...,L/2 obtained with CiSSA.
% I:   Four options:
%      1) A positive integer. It is the number of data per year in
%      economic time series. The function automatically computes the
%      trend (oscillations with period greater than 8 years), the
%      business cycle (oscillations with period between 8 & 1.5 years)
%      and seasonality.
%      2) A cell array. Each cell contains a row vector with the desired
%      values of k to be included in a group, k=1,2,...,L/2. The function
%      computes the reconstructed components for these groups.
%      3) A number between 0 & 1. This number represents the accumulated
%      share of the psd achieved with the sum of the share associated to
%      the largest eigenvalues. The function computes the original
%      reconstructed time series as the sum of these components.
%      4) A number between -1 & 0. It is a percentile (in positive) of
%      the psd. The function computes the original reconstructed time
%      series as the sum of the reconstructed componentes by frequency
%      whose psd is greater that this percentile.
%
% Output arguments:
% rc:  Matrix whose columns are the reconstructed components for each
%      group or subset of frequencies. In the case of economic time series
%      the trend, business cycle and seasonality are in the first, second
%      and third columns, respectively.
% sh:  Column vector with the share(%) of the psd for each group.
% kg:  Cell array where each cell contains a row vector with the values
%      of k belonging to a group. Option 1) produces 3 groups, option 2)
%      gives the goups introduced by the user and options 3) and 4) produce
%      a single group. In option 3), the values of k are sorted according
%      to the share in total psd of their corresponding eigenvalues.
%
% See also: cissa

% -------------------------------------------------------
% Checking the input arguments
% -------------------------------------------------------
% Length and number of the reconstruted series
[T, F] = size(Z);

% Window length
L = length(psd);

% Type and value of input argument #3
if iscell(I)
    opc = 2;
elseif isnumeric(I)
    if (I-floor(I))==0&&I>0
        opc = 1;
    elseif (0<I)&&(I<1)
        opc = 3;
    elseif (-1<I)&&(I<0)
        opc = 4;
    else
        error('***  Input argument #3: Value not valid  ***')
    end
else
    error('***  Input argument #3: Type not valid  ***')
end

switch opc
    case 1
        % Proportionality of L
        if mod(L,I)
            error('***  L is not proportional to the number of data per year  ***');
        end
    case 2
        % Number of groups
        G = length(I);
        if G>F
            error('***  The number of groups is greater than the number of frequencies  ***')
        end    
        % Disjoint groups
        for j=1:G-1
            for m=j+1:G
                if intersect(I{j},I{m})
                    error('***  The groups are not disjoint  ***');
                end
            end
        end   
end

% -------------------------------------------------------
% PSD for frequencies <= 1/2
% -------------------------------------------------------
if mod(L,2)
    pzz = [psd(1); 2*psd(2:F)];
else
    pzz = [psd(1); 2*psd(2:F-1); psd(F)];
end

% -------------------------------------------------------
% Indexes k for each group
% -------------------------------------------------------
switch opc
    case 1
        % Number of groups
        G = 3;
        % Number of data per year
        s = I;
        % Inizialitation of cell array
        kg = cell(G,1);
        % Seasonality
        kg{3} = L*(1:s/2)./s+1;
        % Business cycle
        kg{2} = max(2,floor(L/(8*s)+1)):min(F,floor(L/(1.5*s)+1));
        % Trend
        kg{1} = 1:kg{2}(1)-1;
    case 2
        % Groups
        kg = I;
    case 3
        % Number of groups
        G = 1;
        % Inizialitation of cell array
        kg = cell(G,1);
        % Eigenvalues in decreasing order
        [psor, ks] = sort(pzz,'descend');
        % Cumulative share in percentage
        pcum = 100*cumsum(psor)/sum(psd);
        % Group for the reconstructed time series
        kg{1} = ks(1:length(ks(pcum<100*I))+1);
    case 4
        % Number of groups
        G = 1;
        % Inizialitation of cell array
        kg = cell(G,1);
        % All k values
        ks = 1:F;
        % Group for the reconstructed time series
        kg{1} = ks(pzz>prctile(pzz,-100*I));
end

% -------------------------------------------------------
% Output arguments
% -------------------------------------------------------
% Inizialitation of output arguments
rc = zeros(T,G);
sh = zeros(G,1);

% Computing output arguments
for j=1:G
    % Reconstructed component for each group
    rc(:,j) = sum(Z(:,kg{j}),2);
    
    % Psd share for each group
    sh(j) = 100*sum(pzz(kg{j}))/sum(pzz);
end
