%% RELOAD TRIAL METADATA
% Prasanna Sritharan, October 2019
%
% Reload MAT file


function aclr = loadMetaStruct()

    % get user settings
    addpath('..');
    user = getUserScriptSettings();

    % assign struct fields
    outpath = user.OUTPATH;
    
    % load lpaclr struct
    dbase = load([outpath '\aclr.mat']);
    aclr = dbase.aclr;
        
end