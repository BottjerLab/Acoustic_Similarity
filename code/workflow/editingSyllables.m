function approvedSyllables = editingSyllables(birdID, sessionID)
%Jenny added
% this script is only for editing/adjusting/approving/tweaking
%function tmpSyllables = editingSyllables(birdID, sessionID)
if nargin < 1
    birdID = 'Y231';
end
if nargin < 2
    sessionID = 'Y231_12_12_026';
end
catalog = reportOnData(birdID, sessionID, [], 'verbose',false);

if isempty(catalog), %Jenny changed from ~isempty
    error('Session ID %s does not exist...',sessionID);
end

fprintf('Manifest contents...\n');
{catalog.manifest.name} %#ok<NOPRT>
%% define the specific regions to operate on (like function arguments)
ROIs = [];
while isempty(ROIs)
    varname = input('Which events would you like to load as the visualization ranges (default is ''manualMotifs'')? ','s');
    if isempty(varname)
        varname = 'manualMotifs';
    end
    ROIs = loadFromManifest(catalog.manifest, varname);
    if isempty(ROIs)
        fprintf('Variable %s is not found, try again...\n',varname);
    elseif ~isEvent(ROIs)
        fprintf('Variable %s is not an event, try again...\n',varname);
        ROIs = [];
    end
end

subROIs = [];
while isempty(subROIs)
    varname = input('Which events would you like to load for editing? ', 's');
    subROIs = loadFromManifest(catalog.manifest, varname);
    if isempty(subROIs)
        fprintf('Variable %s is not found, try again...\n',varname);
    elseif ~isEvent(subROIs)
        fprintf('Variable %s is not an event, try again...\n',varname);
        subROIs = [];
    end
end

% load other values
songStruct = loadFromManifest(catalog.manifest, 'songStruct');
noiseMask  = loadFromManifest(catalog.manifest, 'noiseMask');
params = processArgs(defaultParams, birdID, 'fs', 1/songStruct.interval); %Jenny moved this down here

%% resets new results if you want to start over
startMotif = 1;
approvedSyllables = initEvents;
yn = input('Would you like to load approved syllables from the base workspace, y/n [no by default]? ','s');
if strcmpi(yn, 'y'),    
    approvedSyllables = uigetvar('struct');
    startMotif = input('Which motif will you start from (number, please)? ');
end

fprintf('Notice: approvedSyllables dumps to base workspace on abort or matlab error...\n');
%% do actual editing
%%memFlag = false;

% exclude the TUT only bouts - limited time only
ROIs(strcmp({ROIs.type},'TUT')) = [];

ii = startMotif;
while ii <= numel(ROIs) % if you want to continue, change this range
    thisEv = ROIs(ii);
    fprintf('Editing %d/%d...\n', ii, numel(ROIs));
    
    %%try % handle out of memory errors by reducing the number of frequency bands
    revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
        'editSpecType', 'inter', ... 
        'inter.freqBands', linspace(1, 10240, params.inter.NfreqBands), ... % <--- JS: the range of visualization is now 0-10.2kHz 
        'adjustLabels', true, ...
        'dgram.minContrast',3e-12,... %<--- spectrogram contrast is edited HERE (just visualization)
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    %{
    catch err
        if ~strcmp(err.identifier,'MATLAB:nomem')
      
    % try splitting the ROI according to the sub-ROIs 
            aa = subROIs(isWithinEvent(subROIs,thisEv));
            
            splitLens = [aa(2:end).start] - [aa(1:end-1).stop];
            splitTimes = [aa(1:end-1).stop];
            
            [splitLens; splitTimes]
            keyboard
            
            fprintf('Memory overloaded, splitting bout into 2 parts, length %0.2fs and %0.2fs...\n', 1,1);
            % NB: this is BAD practice - don't be like me, Jenny
            % use while loops instead of this pattern
            continue;
   
    else % abort! but clean up first
            cleanupObj = onCleanup(@() putvar(approvedSyllables)); 
            rethrow(err);
end
    end
    %}
    approvedSyllables = [approvedSyllables revisedSylls'];  %#ok<AGROW>
    cleanupObj = onCleanup(@() putvar(approvedSyllables)); 
    %{
 if memFlag % reset the number of bands
        memFlag = false;
        params.inter.NfreqBands = oldBands;
        fprintf('Resetting number of frequency bands to %d...\n', params.inter.NfreqBands);
    end
 %}   
    ii = ii+1;
end

% save interactively - you can cancel without problems
fileToSave = [pwd filesep 'data' filesep birdID filesep 'approvedSyllables-' sessionID '.mat'];
uisave('approvedSyllables', fileToSave);

end