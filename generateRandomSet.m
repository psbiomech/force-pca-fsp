function [means,percentiles] = generateRandomSet(nobs,nvars)

%GENERATERANDOMSET Generate random set for Horn's Parallel Analysis
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% parameters
ndatsets = 1000;  % number of ndatsets
nthpercent = 95;    % desired percentile

% generate correlation matrix for random variables and extract eigenvalues
for nds=1:ndatsets
    evals(:,nds) = eig(corrcoef(randn(nobs,nvars)));
end

% sort the rows of the eigenvalue matrix in descending order
evals = sort(evals,1,'descend');

% sort the eigenvalues for each column.
evals = sort(evals,2);

% extract the nth percentile of the eigenvalues for each column
percentiles = evals(:,round((nthpercent*ndatsets)/100));

% calculate the mean eigenvalues for each position (column).
means = mean(evals,2);


end