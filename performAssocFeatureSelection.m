function [acorrtable, associdx] = performAssocFeatureSelection(final, training, pcainfo)

%PERFORMASSOCFEATURESSELECTION Find PCs associated with primary features
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;
groups = user.GROUPS;

fprintf('Associated features analysis.\n');
fprintf('------------------------------------------------\n'); 

% Data tables
acorrtable = {};
associdx = {};


% Limbs
for b=1:2

    fprintf('\nAssociated feature selection: %s LIMB\n', upper(limbs{b}));

    % Determine correlation between retained PCs and all PCs after Parallel
    % Analysis
    corrassoc = corr(training.(limbs{b}).data, training.(limbs{b}).data(:, final.(limbs{b}).trainidx));
    
    % Determine associated features
    x = 1;
    acorrtable.(limbs{b}) = {};
    for f=1:length(final.(limbs{b}).trainidx)
    
        % Correlation coefficients
        corrassoc(final.(limbs{b}).trainidx(f), f) = NaN;   % ignore self correlations
        iscorr = abs(corrassoc(:, f)) > 0.5;
        associdx.(limbs{b}){x} = find(iscorr)';
        
        fprintf('%s\n',final.(limbs{b}).labels{f});
        
        % If no associated featues, then record empty row
        if isempty(associdx.(limbs{b}){x})
            fprintf('\t----> None\n');
            acorrtable.(limbs{b}){x, 1} = final.(limbs{b}).labels{f};
            acorrtable.(limbs{b}){x, 2} = 'None';
            x = x + 1;
            continue;
        end
        
        % Build output table
        for j = associdx.(limbs{b}){x}
            
            fprintf('\t----> %s\n', training.(limbs{b}).labels{j});
            
            % Descriptives
            fprintf('\t\t\tCalculating descriptives...\n');
            acorrtable.(limbs{b}){x, 1} = final.(limbs{b}).labels{f};
            acorrtable.(limbs{b}){x, 2} = training.(limbs{b}).labels{j};
            acorrtable.(limbs{b}){x, 3} = mean(training.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), j));
            acorrtable.(limbs{b}){x, 4} = std(training.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), j));
            acorrtable.(limbs{b}){x, 5} = mean(training.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), j));
            acorrtable.(limbs{b}){x, 6} = std(training.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), j));
            
            % Correlation coefficient (rho)
            acorrtable.(limbs{b}){x, 7} = corrassoc(j, f);
            
            % T-test (t and p)
            fprintf('\t\t\tPerforming two-sample t-tests...\n');
            [~, pval, ~, stats] = ttest2(training.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), j), training.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), j));
            acorrtable.(limbs{b}){x, 8} = stats.tstat;
            acorrtable.(limbs{b}){x, 9} = pval;
            
            % Effect size (Hedges g)
            fprintf('\t\t\tPerforming effect size using Hedges g...\n');        
            mesvals = mes(training.(limbs{b}).data(pcainfo.(limbs{b}).(['is' groups{1}]), j), training.(limbs{b}).data(~pcainfo.(limbs{b}).(['is' groups{1}]), j), 'hedgesg');
            acorrtable.(limbs{b}){x, 10} = mesvals.hedgesg;
            
            % Increment counter
            x = x + 1;
            
        end

    end

end

% Save results
save(fullfile(outpath,'assocfeatures.mat'),'acorrtable','associdx');


fprintf('------------------------------------------------\n'); 

end

