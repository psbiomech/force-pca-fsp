function trimmed = trimData(original,val0,valEnd,keycol)


%TRIMDATA Trim data between first and last points
%   Prasanna Sritharan, March 2020
%
% Trim data between val0 and valEnd. Data should be arranged in columns.
% The key column should be sorted in increasing order.
%
% Parameters:
%   original: the original untrimmed data vector or matrix.
%   val0: the first value in the key column representing the start of the
%       trimmed data
%   valEnd: the last value in the key column representing the end of the
%       trimmed data
%   keycol: optional integer indicating which column in the data is the key
%       column, data vectors always have key column 1 (default: 1)


% default arguments
if nargin==3, keycol=1; end


% key column data
key = original(:,keycol);

% get indices of first and last data point in key
idx0 = find(key<=val0,1,'last');
idxEnd = find(key>=valEnd,1,'first');

% check index bounds
if isempty(idx0)||(idx0<1), idx0=1; end
if isempty(idxEnd)||(idxEnd>length(key)), idxEnd=length(key); end

% trim data
trimmed = original(idx0:idxEnd,:);



end

