function [pcsexplained, totalvalidpcs, totalvariance, paselected, paquantl] = performParallelAnalysis(pcadata, pcaweights, pcaout, pcainfo)

%PERFORMPARALLELANALYSIS Undertake Horn's Parallel Analysis for FORCe SDP
%   Prasanna Sritharan, February 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;


fprintf('Horn''s Parallel Analysis on PCA scores.\n');
fprintf('------------------------------------------------\n'); 


% Load the 95th percentiles of the eigenvalues of the normally-distributed
% random variables if it exists, otherwise create one and save to file
eigfile = fullfile(outpath, 'randvareig.mat');
if exist(eigfile,'file')
    fprintf('Random variable file exists. Loading...\n');
    for b=1:2
        fprintf('---> %s%s limb\n', upper(limbs{b}(1)), limbs{b}(2:end));
        aux = load(eigfile);
        randvareigpcts.(limbs{b}) = aux.randvareigpcts.(limbs{b});
    end
else
    fprintf('Random variable file does not exist. Generating random variables...\n');
    for b=1:2
        fprintf('---> %s%s limb\n', upper(limbs{b}(1)), limbs{b}(2:end));
        obs = pcainfo.(limbs{b}).observations.total;
        vars = pcainfo.(limbs{b}).variables;
        [randvareigmeans.(limbs{b}), randvareigpcts.(limbs{b})] = generateRandomSet(obs,vars);
    end
    save(eigfile, 'randvareigmeans', 'randvareigpcts');
end

% Data tables
pcsexplained = struct;
totalvariance = struct;
totalvalidpcs = struct;
paselected = struct;
paquantl = struct;

% Limbs
for b=1:2

    fprintf('\nHorn''s Parallel Analysis: %s LIMB\n', upper(limbs{b}));

    % Perform Parallel Analysis per variable
    fprintf('---> Performing Parallel Analysis...\n');
    inc = 1;

    % Analysis
    pcsexplained.(limbs{b}) = {};
    totalvariance.(limbs{b}) = {};
    totalvalidpcs.(limbs{b}) = 0;
    for d={'ik','id'}
        
        % Retain eigenvalues of R, the weighted correlation matrix of data X,
        % that exceed the equivalent eigenvalue from the random variable set at
        % the 95th percentile
        dataset = pcadata.(limbs{b}).(d{1}); 
        label = pcainfo.(limbs{b}).(d{1}).label;
        for v=1:size(dataset,3)
            
    
            % Find scores
            data = squeeze(dataset(:, :, v));        
            latent = sort(eig(weightedcorrs(data,pcaweights.(limbs{b}))), 'descend');
            numvalidscores = find(latent>randvareigpcts.(limbs{b}), 1, 'last');
                    
            % Record table of results for Parallel Analysis
            % (how much each selected PC explains the associated waveform)
            explained = pcaout.(limbs{b}).(d{1}).explained(:,v);
            varname = pcainfo.(limbs{b}).(d{1}).varnames{v};
            pcsexplained.(limbs{b}){inc,1} = [label '_' varname];
            for bb=1:numvalidscores
                pcsexplained.(limbs{b}){inc,1+bb} = explained(bb)./100;
            end
            totalvariance.(limbs{b}){inc} = sum(cell2mat(pcsexplained.(limbs{b})(inc, 2:end)));
            totalvalidpcs.(limbs{b}) = totalvalidpcs.(limbs{b}) + numvalidscores;
            inc = inc + 1;
            
            % Record the eigenvector matrix (coefficients), and also the scores
            % for the selected PCs for each variable
            paselected.(limbs{b}).(d{1}).(varname).coeff = squeeze(pcaout.(limbs{b}).(d{1}).coeff(:, :, v));
            paselected.(limbs{b}).(d{1}).(varname).score = squeeze(pcaout.(limbs{b}).(d{1}).score(:, 1:numvalidscores, v));        
                    
            % Get quantiles and indices of scores associated with top and
            % bottom quantiles
            paquantl.(limbs{b}).(d{1}).(varname).quantl = quantile(paselected.(limbs{b}).(d{1}).(varname).score, [.025 .25 .50 .75 .975]);
            paquantl.(limbs{b}).(d{1}).(varname).bottom.idx = paselected.(limbs{b}).(d{1}).(varname).score <= paquantl.(limbs{b}).(d{1}).(varname).quantl(2, :);
            paquantl.(limbs{b}).(d{1}).(varname).top.idx = paselected.(limbs{b}).(d{1}).(varname).score >= paquantl.(limbs{b}).(d{1}).(varname).quantl(4, :);        
            
            % Get the data associated with the PC quantiles
            for bb=1:numvalidscores
                paquantl.(limbs{b}).(d{1}).(varname).bottom.mean(:, bb) = mean(pcadata.(limbs{b}).(d{1})(paquantl.(limbs{b}).(d{1}).(varname).bottom.idx(:,bb), :, v), 1);
                paquantl.(limbs{b}).(d{1}).(varname).bottom.std(:, bb) = std(pcadata.(limbs{b}).(d{1})(paquantl.(limbs{b}).(d{1}).(varname).bottom.idx(:,bb), :, v), [], 1);        
                paquantl.(limbs{b}).(d{1}).(varname).top.mean(:, bb) = mean(pcadata.(limbs{b}).(d{1})(paquantl.(limbs{b}).(d{1}).(varname).top.idx(:,bb), :, v), 1);
                paquantl.(limbs{b}).(d{1}).(varname).top.std(:, bb) = std(pcadata.(limbs{b}).(d{1})(paquantl.(limbs{b}).(d{1}).(varname).top.idx(:,bb), :, v), [], 1);
                paquantl.(limbs{b}).(d{1}).(varname).normpc(:, :, bb) = weightedcorrs([paselected.(limbs{b}).(d{1}).(varname).score(:,bb) squeeze(pcadata.(limbs{b}).(d{1})(:, :, v))], pcaweights.(limbs{b}));
                paquantl.(limbs{b}).(d{1}).(varname).normpc2(:, bb) = squeeze(paquantl.(limbs{b}).(d{1}).(varname).normpc(2:end, 1, bb)).^2;
            end
            
        end   
    
    end

end

% Save results
fprintf('\nSaving Parallel Analysis inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'parallelanalysis.mat'),'pcsexplained','totalvalidpcs','totalvariance','paselected','paquantl');

fprintf('------------------------------------------------\n'); 

end

