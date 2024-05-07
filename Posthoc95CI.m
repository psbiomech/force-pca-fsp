%POSTHOC95CI Calculate 95% CI of mean differences for T-table
%   Prasanna Sritharan, April 2024
%
% This should be undertaken as part of the t-table generation, but has been
% added post-hoc in response to reviewer request during peer review of the
% RSOS/RSIF paper associated with this work.


% ## DATA FROM TABLES

% Load data
srcpath = 'C:\Users\Owner\Documents\data\FORCe\pcaDatabase\by-leg\tables';
data_pivot = readcell(fullfile(srcpath, 'table2_ttestout_pivot.csv'));
data_nonpivot = readcell(fullfile(srcpath, 'table2_ttestout_nonpivot.csv'));
data = {data_pivot, data_nonpivot};

% Sample sizes (no of trials, symptomatic vs control)
ns = [193 197];

 
% Calculate 95% CI of mean difference
leg = {'pivot', 'nonpivot'};
for d=1:2

    % Storage
    md95 = cell(11,1);
    ci95_upper = cell(11,1);
    ci95_lower = cell(11,1);

    % Column headers
    md95{1} = 'mean_diff';
    ci95_upper{1} = '95pc_ci_upper';
    ci95_lower{1} = '95pc_ci_lower';

    % Evaluate each feature
    for f=1:10

        % Data row
        r = f + 1;

        % Group descriptives
        ms = [data{d}{r,2}, data{d}{r,4}];
        sds = [data{d}{r,3}, data{d}{r,5}];

        % Mean difference
        md = ms(1) - ms(2);

        % Pooled standard deviation
        dfs = ns(1) + ns(2) - 2;
        psd = sqrt(((ns(1)-1)*(sds(1)^2) + (ns(2)-1)*(sds(2)^2)) / dfs);

        % T-score (two-sided p=0.05)
        tscore = tinv(0.975, dfs);

        % CI range
        cirange = tscore * psd * sqrt( (1/ns(1)) + (1/(ns(2))));

        % 95% CI of the mean difference
        ci95md = [md-cirange, md+cirange];
   
        % Store
        md95{r} = md;
        ci95_lower{r} = ci95md(1);
        ci95_upper{r} = ci95md(2);

    end

    % update table
    data_updated = [data{1}, md95, ci95_lower, ci95_upper];

    % Save results
    fprintf(['Saving T-tests inputs and outputs (', leg{d},')...\n']);
    if ~exist(srcpath,'dir'), mkdir(srcpath); end
    writecell(data_updated, fullfile(srcpath, ['ttable_with_95ci_', leg{d}, '.csv']));

end
