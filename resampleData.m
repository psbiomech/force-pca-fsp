%% RESAMPLE ONTARIO DATA
% Prasanna Sritharan, October 2012
%
% Resample a vector or matrix. For matrices, data should be arranged in
% columns
%
% invector: Data arranged in columns. Each row is a time step.
% RESAMP: desired number of samples




function outvector = resampleData(invector,RESAMP)

% output vector
outvector0 = zeros(RESAMP,size(invector,2));

% length of input data vector
inlength = size(invector,1);

% pad vector to prevent edge effects
paddedin = [repmat(invector(1,:),inlength,1); invector; repmat(invector(end,:),inlength,1)];

% resample column by column
for i=1:size(invector,2)

    % resample padded vector
    paddedout = resample(paddedin(:,i),3*RESAMP,3*inlength);

    % output resampled vector
    outvector0(:,i) = paddedout(RESAMP+1:2*RESAMP);    

end

% return resampled data
outvector = outvector0;

end