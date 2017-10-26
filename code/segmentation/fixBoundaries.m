function newRegions = fixBoundaries(songStruct, regions, params, varargin)

newRegions = initEvents(numel(regions));
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});
noRollParams = processArgs(params, 'preroll',0,'postroll',0);
% depends heavily on params.syllable.minPower as that is the final
% threshold (in band pass)
fs = 1/songStruct.interval;
margin = params.boundaryAdvance; % seconds, since this is a slow iterative algorithm, this is a very performance-affecting parameter
silentMinLength = params.minimumSilentGap;
incPrerollParams = processArgs(params,'preroll', 0, 'postroll',0); % milliseconds
incPostrollParams = processArgs(params,'postroll', 0, 'preroll',0); % milliseconds

for ii = 1:numel(regions)
    
    fprintf('Repairing region %d...\n',ii);
    thisRegion = regions(ii);

    % start
    regionPriorStart = thisRegion.start;
    newStart = regionPriorStart;
    onsetBracket = eventFromTimes(thisRegion.start, thisRegion.start + 5 * margin, fs);
    while newStart == regionPriorStart
        onsetBracket = addPrePost(onsetBracket, incPrerollParams);
        silentRegions = findSilence(songStruct, onsetBracket, noRollParams);

        % make the preroll a little bigger each time
        incPrerollParams.preroll = incPrerollParams.preroll + 1000 * margin;
            
        % if there is longer than a margin duration of unbroken silence
        % that begins BEFORE or AT the prior region start
        silBeforeStart = [silentRegions.start] <= regionPriorStart;
        silSuffLength = ([silentRegions.stop] - [silentRegions.start]) > silentMinLength;
        if any(silBeforeStart & silSuffLength)
            latestReg = silentRegions(...
                find(silBeforeStart & silSuffLength,1,'last'));
            newStart = latestReg.stop - margin;
        end
    end
    

    % stop
    regionPriorStop = thisRegion.stop;
    newStop = regionPriorStop;
    offsetBracket = eventFromTimes(thisRegion.stop - 5 * margin, thisRegion.stop, fs);
    while newStop == regionPriorStop
        offsetBracket = addPrePost(offsetBracket, incPostrollParams);
        silentRegions = findSilence(songStruct, offsetBracket, noRollParams);

        % make the preroll a little bigger each time
        incPostrollParams.postroll = incPostrollParams.postroll + 1000 * margin;
            
        % if there is longer than a margin duration of unbroken silence
        % that ends AFTER or AT the prior region start
        silAfterStop = [silentRegions.stop] >= regionPriorStop;
        silSuffLength = ([silentRegions.stop] - [silentRegions.start]) > silentMinLength;
        if any(silAfterStop & silSuffLength)
            earliestReg = silentRegions(...
                find(silAfterStop & silSuffLength,1,'first'));
            newStop = earliestReg.start + margin;
        end
    end
      
    newRegions(ii) = eventFromTimes(newStart, newStop, fs);
% test only
    if params.plot
        [silentRegions, soundRegions, spec] = findSilence(songStruct, newRegions(ii), noRollParams);

        plotAllFigures(spec, soundRegions, params);
        plotMarks(adjustTimeStamps(thisRegion,-newRegions(ii).start));
    end
    if params.verbose
          fprintf('Adding %.2fs to the onset, %.2fs to the offset...\n', ...
                    thisRegion.start - newRegions(ii).start , ...
                    newRegions(ii).stop - thisRegion.stop);
    end
end

