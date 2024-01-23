 function user = getUserScriptSettings()

%GETUSERSCRIPTSETTINGS User settings for processing
%   Prasanna Sritharan, June 2022



%% FOLDER PATHS


% folder paths: Lenovo laptop
user.CODEROOT = 'C:\Users\Owner\Documents\data\FORCe\'; 
user.OUTPATH1 = 'C:\Users\Owner\Documents\data\FORCe\pcaDatabase\by-leg\';   % location of output database
user.SRCPATH = 'C:\Users\Owner\Documents\data\FORCe\outputDatabase\csvfolder';   % location of source data



%% GENERAL
% ------------------------------

% general parameters
user.LIMBS = {'pivot', 'nonpivot'};
user.GROUPS = {'sym','ctrl'};
user.SUBJPREFIX = {'FAILT', 'FAILTCRT'};
user.TRIALCOMBO = {{'pivot_more', 'pivot'}, {'pivot_less', 'pivot'}};
user.FOOT = {'r','l'};
user.RESAMP = 101;
user.GRAVITY = 9.81;    % m/s2
user.DATACOLS = 32:132;



%% FEATURE SELECTION PARAMETERS: BY LEGS

% Rajagopal model
user.feature.rajagopal.ik.label = 'angle';
user.feature.rajagopal.ik.headers = {'hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','lumbar_extension','lumbar_bending','lumbar_rotation'};
user.feature.rajagopal.ik.idx = [8:10 11 13 16:18];

user.feature.rajagopal.id.label = 'moment';
user.feature.rajagopal.id.headers = {'hip_flexion','hip_adduction','hip_rotation_moment','knee_angle_moment','ankle_angle_moment','lumbar_extension_moment','lumbar_bending_moment','lumbar_rotation_moment'};
user.feature.rajagopal.id.idx = [8:10 14 19 11:13];


%% FEATURE SELECTION PARAMETERS: BOTH LEGS

% Rajagopal model
user.feature2.rajagopal.ik.label = 'angle';
user.feature2.rajagopal.ik.headers = {'hip_adduction','hip_flexion','hip_rotation','knee_angle','ankle_angle','lumbar_bending','lumbar_extension','lumbar_rotation', 'hip_adduction_np','hip_flexion_np','hip_rotation_np','knee_angle_np','ankle_angle_np'};
user.feature2.rajagopal.ik.idx = [8:10 11 13 16:18 33:35 36 38];

user.feature2.rajagopal.id.label = 'moment';
user.feature2.rajagopal.id.headers = {'hip_adduction','hip_flexion','hip_rotation','knee_angle','ankle_angle','lumbar_bending','lumbar_extension','lumbar_rotation', 'hip_adduction_np','hip_flexion_np','hip_rotation_np','knee_angle_np','ankle_angle_np'};
user.feature2.rajagopal.id.idx = [8:10 14 19 11:13 33:35 39 44];


end
