function [wavecorr, dataidx] = performWaveformCorr(final, pcainfo, pcadata, pcaweights)


%PERFORMPCWAVEFORMCORR Correlation of PC scores against original waveforms
%   Prasanna Sritharan, March 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;

fprintf('Correlations of PC scores against original waveform data.\n');
fprintf('------------------------------------------------\n'); 

% Data tables
wavecorr = struct;
dataidx = struct;

% Limbs
labels = {'ik','id'};
alias = {'angle','moment'};
pcregexp = '(moment|angle)_(\w+)_PC\d+$';
for b=1:2

    fprintf('\nWaveform correlations: %s LIMB\n', upper(limbs{b}));

    % Perform correlation against waveforms
    dataidx.(limbs{b}) = cell(length(final.(limbs{b}).labels), 2);
    for s=1:length(final.(limbs{b}).labels)
        
        fprintf('%s\n', final.(limbs{b}).labels{s});
        
        % Get indices of original waveform data
        tokens = regexpi(final.(limbs{b}).labels{s}, pcregexp, 'tokens');
        dataidx.(limbs{b}){s, 1} = labels{strcmpi(alias, tokens{1}{1})};
        dataidx.(limbs{b}){s, 2} = find(strcmpi(pcainfo.(limbs{b}).(dataidx.(limbs{b}){s, 1}).varnames, tokens{1}{2}), 1);
    
        % Calculate waveform correlation coefficients
        fprintf('---> Calculating correlation coefficients...\n');
        wavecorr.(limbs{b}).coeffs(:, :, s) = weightedcorrs([final.(limbs{b}).data(:,s) squeeze(pcadata.(limbs{b}).(dataidx.(limbs{b}){s, 1})(:, :, dataidx.(limbs{b}){s, 2}))], pcaweights.(limbs{b}));
        wavecorr.(limbs{b}).norms(:, s) = squeeze(wavecorr.(limbs{b}).coeffs(2:end,1,s)).^2; 
        
    end

end

% Save results
fprintf('\nSaving Waveform Correlation inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'waveformcorr.mat'),'wavecorr','dataidx');


fprintf('------------------------------------------------\n'); 

end

