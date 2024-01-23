function tbls = generateOutputTables(pcsexplained, totalvariance, ttable, acorrtable)


%GENERATEOUTPUTTABLES Generate output tables
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;

fprintf('Generating output tables.\n');
fprintf('------------------------------------------------\n');

% Data tables
tbls = {};

% Limbs
for b = 1:2

    tbls.(limbs{b}) = {};

    fprintf('\nTables: %s LIMB\n', upper(limbs{b}));

    % Output folder
    tabledir = fullfile(outpath, 'tables');
    if ~exist(tabledir, 'dir'), mkdir(tabledir); end
    
    % Table 1: variance explained by PCs
    fprintf('Creating Table 1: Variance explained by PCs...\n');
    pcsexplained_percent_1decpl = pcsexplained.(limbs{b});
    pcsexplained_percent_1decpl(:, 2:end) = cellfun(@(x) round(x*100, 1), pcsexplained_percent_1decpl(:, 2:end), 'UniformOutput', false);
    totalvariance_percent_1decpl = cellfun(@(x) round(x*100, 1), totalvariance.(limbs{b}), 'UniformOutput', false);
    headers1 = ['variable', cellfun(@(x) ['PC' num2str(x)], num2cell(1:size(pcsexplained_percent_1decpl, 2) - 1), 'UniformOutput', false), 'total'];
    tbls.(limbs{b}).table1 = [headers1; [pcsexplained_percent_1decpl, totalvariance_percent_1decpl']];
    writecell(tbls.(limbs{b}).table1, fullfile(tabledir, ['table1_pcsexplained_' limbs{b} '.csv']));
    
    % Table 2: t-test, effect size and descriptives
    fprintf('Creating Table 2: T-test, effect size and descriptives for final selected PCs...\n');
    headers2 = {'feature', 'mean_sym', 'sd_sym', 'mean_ctrl', 'sd_ctrl', 't', 'P', 'g'};
    tbls.(limbs{b}).table2 = [headers2; ttable.(limbs{b}).table];
    writecell(tbls.(limbs{b}).table2, fullfile(tabledir, ['table2_ttestout_' limbs{b} '.csv']));
    
    % Table 3: qualitative interpretation of PCs
    % N/A, descriptive table
    fprintf('Creating Table 3: Qualitative interpretation of PCs...\n');
    tbls.(limbs{b}).table3 = 'N/A, qualitative interpretation of PCs, descriptive table';
    
    % Table 4: associated features
    fprintf('Creating Table 4: Associated features...\n');
    headers4 = {'main_feature', 'associated_feature', 'mean_sym', 'sd_sym', 'mean_ctrl', 'sd_ctrl', 'rho', 't', 'P', 'g'};
    tbls.(limbs{b}).table4 = [headers4; acorrtable.(limbs{b})];
    writecell(tbls.(limbs{b}).table4, fullfile(tabledir, ['table4_assocfeatures_' limbs{b} '.csv']));

end


% Save results
fprintf('\nSaving output tables...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'output.mat'),'tbls');


fprintf('------------------------------------------------\n');

end

