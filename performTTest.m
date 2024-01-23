function ttable = performTTest(final, pcainfo)


%PERFORMTTEST Perform T-test on final feature set
%   Prasanna Sritharan, March 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;
groups = user.GROUPS;


fprintf('T-tests, effect size and descriptives for final PC set.\n');
fprintf('------------------------------------------------'); 

% Data tables
ttable = struct;

% Limbs
for b=1:2

    fprintf('\nT-Tests: %s LIMB\n', upper(limbs{b}));

    % Run tests and build output table
    npcs = length(final.(limbs{b}).labels);
    ttable.(limbs{b}).table = cell(npcs, 8);
    for n=1:npcs
        
        fprintf('%s\n',final.(limbs{b}).labels{n});
        
        % Descriptives
        fprintf('---> Calculating descriptives...\n');
        ttable.(limbs{b}).table{n,1} = final.(limbs{b}).labels{n};
        ttable.(limbs{b}).table{n,2} = mean(final.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), n));
        ttable.(limbs{b}).table{n,3} = std(final.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), n));
        ttable.(limbs{b}).table{n,4} = mean(final.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), n));
        ttable.(limbs{b}).table{n,5} = std(final.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), n));    
        
        % T-test (t and p)
        fprintf('---> Performing two-sample t-tests...\n');
        [~, pval, ~, stats] = ttest2(final.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), n), final.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), n));
        ttable.(limbs{b}).table{n,6} = stats.tstat;
        ttable.(limbs{b}).table{n,7} = pval;
        
        % Effect size (Hedges g)
        fprintf('---> Performing effect size using Hedges g...\n');
        mesout = mes(final.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), n), final.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), n), 'hedgesg');
        ttable.(limbs{b}).table{n,8} = mesout.hedgesg;
        
        % Headers
        %ttable.(limbs{b}).headers = {'feature', 'mean1', 'sd1', 'mean2', 'sd2', 't', 'p', 'g'};

    end
    
end

% Save results
fprintf('\nSaving T-tests inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'ttable.mat'),'ttable');

fprintf('------------------------------------------------\n');

end

