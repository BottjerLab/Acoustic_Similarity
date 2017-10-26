function editedRegions = plotAndAdjust(songStruct, changeRegions, parentRegion, params, varargin)
%PLOTANDADJUST Plots for manual event adjustment
%
%   function editedRegions = plotAndAdjust(songStruct, changeRegions, regionRange)
%   acquires metrics about the audio clip in regionRange(s), and plots the 
%   information in a console.  It plots the regions in changeRegions and  
%   allows user control over the region boundaries, as described in editEventsLabel.  
%
%   Parameter options:
%   optGraphs: controls what information will be plotted
%   adjustLabels: determines if labels will be adjusted
%   editSpecType: determines the specificity of the spectrogram 
%   (look in defaultParams.m for possible options)

if nargin < 4 || isempty(params)
    params = defaultParams;
end

params = processArgs(params,varargin{:});

% batch mode for multiple region ranges
if numel(parentRegion) > 1 
    disp('Multiple regions specified, running on one at a time:');
    editedRegions = initEvents;
    for ii = 1:numel(parentRegion)
        editedRegions = [editedRegions plotAndAdjust(songStruct, changeRegions, parentRegion(ii), params)];
    end
    return;
end

% now regionRange should always be a single entry, check if it is a region
if ~isEvent(changeRegions) || ~isEvent(parentRegion)
    error('plotAndAdjust:badArgs','Arguments should be regions');
end
    
fs = 1/songStruct.interval;
params.(params.editSpecType).fs = fs;
params.fs = fs;

% get the sub regions
changeRegions = getSubEvents(parentRegion, changeRegions);
% verify the regions
[isValid, repairedRegions] = checkRegions(changeRegions, fs);
if any(~isValid)
    ch = input('Region indices and times are inconsistent, repair (Y/n)? ', 's');
    if ~isempty(ch), ch = lower(ch(1)); end;
    if strcmp(ch,'n')
        warning('Regions not repaired, may yield wrong results...');
    else
%         ch = input('Manual repair (y/N)? ', 's');
%         if ~isempty(ch), ch = lower(ch(1)); end;
%         if strcmp(ch,'y')
%             keyboard
%         else
            changeRegions = repairedRegions;
%         end
    end
end
editedRegions = changeRegions;

% get the clip with processing
wholeRangeCl = getClipAndProcess(songStruct, parentRegion, params);
wholeRange = addPrePost(parentRegion, params);

% get all the statistics about a clip
if ~all(strcmp(params.optGraphs,'waveform'))
    wholeRangeSpec = getMTSpectrumStats(wholeRangeCl, params.(params.editSpecType));
else
    wholeRangeSpec.waveform = wholeRangeCl;
    wholeRangeSpec.times = [0 numel(wholeRangeCl)/fs];
end

% get the regions
if ~isempty(wholeRangeSpec)
    % plot the region to be changed within the entire clip
    changeRegions = adjustTimeStamps(changeRegions, -wholeRange.start,fs);
    plotAllFigures(wholeRangeSpec,changeRegions,params);
    % get the edited version if we're meant to edit
    if ~isfield(params, 'readOnly') || ~params.readOnly
        editedRegions = editEventsLabel(changeRegions,fs,params.adjustLabels);
        editedRegions = adjustTimeStamps(editedRegions, +wholeRange.start,fs);
    end
else
    editedRegions = changeRegions;
end

% adopt parent label
if ~isempty(editedRegions) && ~params.adjustLabels
    [editedRegions.type] = deal(parentRegion.type);
end
end