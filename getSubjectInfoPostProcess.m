function subjlist = getSubjectInfoPostProcess(aclr)


%GETSUBJECTINFOPOSTPROCESS Demographic data for publication
%   Prasanna Sritharan, Mario Andrés Muñoz, April 2021
%
% For post-processing only. The function performWeightedPCA() has been
% updtaed to include this functionality.


% user settings
user = getUserScriptSettings();
outpath = user.FEATPATH;
groups = user.GROUPS;
dynamicprefix = user.DYNAMICPREFIX;


% groups
fprintf('Extracting OpenSim data from files...\n');
qs = [];
x = 1;
subjlist = {{},{}};
for g=1:length(groups)

    % subjects
    subjs = fieldnames(aclr.(groups{g}));
    for s=1:length(subjs)
                 
        % trials
        q = 0;
        trials = fieldnames(aclr.(groups{g}).(subjs{s}));
        for t=1:length(trials)        

            % check if struct field is really a trial
            if ~startsWith(trials{t},dynamicprefix), continue; end         

            % ignore trial if flag set
            if aclr.(groups{g}).(subjs{s}).(trials{t}).ignore, continue; end   
            
            % trial folder location and file name prefix
            dest = aclr.(groups{g}).(subjs{s}).(trials{t}).fullpath; 
            fprefix = aclr.(groups{g}).(subjs{s}).(trials{t}).fileprefix;
            
            % trial leg
            leg = aclr.(groups{g}).(subjs{s}).(trials{t}).leg;
            
            % open trial MAT file
            tdata = load(fullfile(dest,fprefix));           
            
            % increment trial counters
            q = q + 1;
            
            % record subject name
            subjlist{g} = [subjlist{g}, subjs{s}];
            
        end
                   
    end
        
end


end

