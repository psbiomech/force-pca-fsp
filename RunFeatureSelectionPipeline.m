%% FEATURE SELECTION PIPELINE - FORCE SDP
% Prasanna Sritharan, May 2022
%
% Note: back up the results of any previous run of this script as the
% Sequential Feature Selection analysis involves random selection of
% principal components, and therefore results from the new run can
% differ substantially from those obtained on a previous run.

% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;

% improve pseudo-randomness by setting seed to time
rng('shuffle');


%% FEATURE SELECTION PIPELINE

% Weighted PCA
[pcadata, pcaweights, pcaout, pcainfo] = performWeightedPCA('rajagopal');

% Parallel Analysis
[pcsexplained, totalvalidpcs, totalvariance, paselected, paquantl] = performParallelAnalysis(pcadata, pcaweights, pcaout, pcainfo);

% Feature selection
[final, weissind, training] = performFeatureSelection(paselected, pcainfo);

% t-tests, effect size and descriptives
ttable = performTTest(final, pcainfo);

% correlation of PCs with original waveforms
[wavecorr, dataidx] = performWaveformCorr(final, pcainfo, pcadata, pcaweights);

% associated features
[acorrtable,associdx] = performAssocFeatureSelection(final,training,pcainfo);


%% TABLES AND FIGURES

% generate output tables
tbls = generateOutputTables(pcsexplained, totalvariance, ttable, acorrtable);

% generate output figures
generateOutputFigures(paselected, paquantl, final);


%% SAVE OUTPUTS

pcaweighted.pcadata = pcadata;
pcaweighted.pcaweights = pcaweights;
pcaweighted.pcaout = pcaout;
pcaweighted.pcainfo = pcainfo;

parallelanalysis.pcsexplained = pcsexplained;
parallelanalysis.totalvalidpcs = totalvalidpcs;
parallelanalysis.totalvariance = totalvariance;
parallelanalysis.paselected = paselected;
parallelanalysis.paquantl = paquantl;

featureselection.final = final;
featureselection.weissind = weissind;
featureselection.training = training;

ttestsig.ttable = ttable;

waveformcorr.wavecorr = wavecorr;
waveformcorr.dataidx = dataidx;

assocfeat.acorrtable = acorrtable;
assocfeat.associdx = associdx;

output.tables = tbls;

% save all structs
save(fullfile(outpath, 'fsp_all_outputs.mat'), 'pcaweighted', 'parallelanalysis', 'featureselection', 'ttestsig', 'waveformcorr', 'assocfeat', 'output');
