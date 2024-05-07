function [final, weissind, training] = performFeatureSelection(paselected, pcainfo)



%PERFORMFEATURESELECTION Perform Feature Selection for FORCe SDP
%   Prasanna Sritharan, February 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)
%
% 1. Build training dataset
% 2. Get upper and lower quantiles for later interpretation of PCs
% 2. Reduce using Weiss-Indurkhya independent feature selection
% 3. Further reduce using Sequential Feature Selection using Naive Bayes
%       Classifier


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;
groups = user.GROUPS;


fprintf('Perform feature selection on IK and ID data.\n');
fprintf('------------------------------------------------'); 

% Data tables
training = struct;
weissind = struct;
final = struct;

% Limbs
for b=1:2

    fprintf('\nFeature selection: %s LIMB\n', upper(limbs{b}));

    % ********************
    % Training data matrix
    
    % prepare training data matrix from selected PCs
    fprintf('---> Preparing training data matrix from PCs selected using Parallel Analysis...\n');
    training.(limbs{b}).data = [];
    training.(limbs{b}).labels = {};
    for d={'ik','id'} 
        label = pcainfo.(limbs{b}).(d{1}).label;
        for v=1:length(pcainfo.(limbs{b}).(d{1}).varnames)        
            varname = pcainfo.(limbs{b}).(d{1}).varnames{v};        
            nsel = size(paselected.(limbs{b}).(d{1}).(varname).score,2);
            training.(limbs{b}).data = [training.(limbs{b}).data, paselected.(limbs{b}).(d{1}).(varname).score];
            training.(limbs{b}).labels = [training.(limbs{b}).labels, cellfun(@(x) [label '_' varname '_PC' num2str(x)], num2cell(1:nsel), 'UniformOutput', false)];    
        end    
    end
    
    
    % ********************
    % Weiss-Indurkhya independent feature selection method
    
    % Get signficance values from Weiss-Indurkhya method
    fprintf('---> Reducing the number of PCs using the Weiss-Indurkhya independent feature selection method...\n');
    indfeatout = IndFeat(training.(limbs{b}).data, pcainfo.(limbs{b}).(['is' groups{1}]));
    
    % Retain only significant PCs
    issigpc = indfeatout>=2.0;  % t-scores
    weissind.(limbs{b}).data = training.(limbs{b}).data(:, issigpc);
    weissind.(limbs{b}).labels = training.(limbs{b}).labels(issigpc);
    weissind.(limbs{b}).trainidx = find(issigpc);
    
    
    % ********************
    % Sequential Feature Selection using Naive Bayes Classifier
    
    fprintf('---> Determining final set of PCs using Sequential Feature Selection...\n');
    
    % Function: number of misclassifications by Naive Bayes Classifier model
    fun = @(XT,yT,Xt,yt)(sum(~eq(yt, predict(fitcnb(XT, yT), Xt))));
    
    % Perform Sequential Feature Selection
    niter = 1000;
    ncSFS = zeros(1, size(weissind.(limbs{b}).data,2));
    opts = statset('display', 'iter', 'UseParallel', false);    % set parallel processing off for local machine as parfor is slower than for loops
    fprintf('\n*** Number of iterations: %d ***\n',niter);
    for n=1:niter
        fprintf('\nITERATION: %d\n',n);
        cvpart = cvpartition(pcainfo.(limbs{b}).(['is' groups{1}]), 'k', 10);    % create k new partitions on every iteration
        fs = sequentialfs(fun, weissind.(limbs{b}).data, pcainfo.(limbs{b}).(['is' groups{1}]), 'cv', cvpart, 'NFeatures', 10, 'Direction', 'forward', 'options', opts);
        ncSFS(fs) = ncSFS(fs) + 1;    % increment counter
    end
    
    % Retain only selected PCs
    fprintf('---> Retaining selected PCs only...')
    [~, pcidx] = sort(ncSFS, 'descend');
    pcidx = sort(pcidx(10:-1:1));
    final.(limbs{b}).data = weissind.(limbs{b}).data(:, pcidx);
    final.(limbs{b}).labels = weissind.(limbs{b}).labels(pcidx);
    final.(limbs{b}).trainidx = weissind.(limbs{b}).trainidx(pcidx);

end

% Save results
fprintf('\nSaving Feature Selection inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'featureselection.mat'),'training','weissind','final');

fprintf('------------------------------------------------\n');

end

