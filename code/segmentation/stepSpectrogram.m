function regions = stepSpectrogram(songStruct, params, varargin)
%STEPSPECTROGRAM Segments a chronic recording to find candidate sounds
%
%  function regions = stepSpectrogram(songStruct) takes in the structure 
%  exported by Spike2, which is expected to have a values field and an
%  interval field defining the sampling rate (as in the interval between 
%  samples).  The function divides the clip into many overlapping subclips
%  and performs a simple segmentation based on volume to find audio regions
%  of interest.  Regions is an array of event structures, with defined
%  onsets/offsets of sounds in both time and sample coordinates.  
%
%   Parameter options:
%   Nsplits: controls what information will be plotted
%   overlapSplit: determines if labels will be adjusted
%   rough: determines the specificity of the spectrogram (look in
%   features/getSpectrumStats for more details)
%    
%   quiet: [t/f] controls whether whole windows are played
%   playsample: [t/f] controls whether segmented regions are played
%   plot: [t/f] controls plotting output
%   verbose: [t/f] controls printing output to command window
%   (look in defaultParams.m for possible options)


if ~exist('params') || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});
params.fs = 1 / songStruct.interval;
fs = params.fs;
regions = initEvents;

% get the clip segment boundaries
wholeSong = getWholeSongAsRegion(songStruct);
windows = splitIntoOverlap(wholeSong, params.Nsplits, params.overlapSplit / 1000);

if ~params.plot
    progressbar('Running rough spectrograms');
end
for ii = 1:params.Nsplits
    if params.verbose
        fprintf('Running segment %d (%s--%s)...\n', ii, ...
            timeToString(windows(ii).start), timeToString(windows(ii).stop));
    end
    sample = getClip(windows(ii), songStruct);
    
    % filter sample with high pass @ 500 Hz   
    params.coarse.fs = 1/songStruct.interval;
    sample = highPassSample(sample,params);

    % play the music
    if ~params.quiet, player = playSound(sample, params.fs); end

    % acquire spectrogram from MATLAB, frequencies, times, psd, characteristic,
    % intercepts, center frequency, and total power
    spec = getMTSpectrumStats(sample, params.coarse);

    % find regions with more than just noise, using simple thresholding
    newRegions = findPossibleSounds(spec, params);
    
    if ~all(checkRegions(newRegions,fs))
        disp('bad newRegions');
        keyboard
    end
    % plotting for visualization
    if params.plot, 
        [hax, optGs, hfig] = plotAllFigures(spec, newRegions, params); 
        title(sprintf('PSD from %s--%s with centroid freq.', ...
            timeToString(windows(ii).start), ...
            timeToString(windows(ii).stop )));

        % draw the appropriate thresholds
        powerPlotHandle = hax(strcmp('totalPower', optGs));
        if ~isempty(powerPlotHandle)
            hold(powerPlotHandle, 'on');
            xx = xlim(powerPlotHandle);
            semilogy(powerPlotHandle, xx, params.powerThresh * [1 1],'b--');
            semilogy(powerPlotHandle, xx, params.riseThresh  * [1 1],'y--');
            title(powerPlotHandle, [get(get(gca,'Title'),'String'), ...
                '  {\color{blue}power \color{yellow}rise}']);
            hold(powerPlotHandle, 'off');
        end

        change_current_figure(hfig);
        drawnow;
    else
        progressbar(ii/params.Nsplits);
    end
    
    % correct regions for true time (in the frame of the whole songStruct)
    newRegions = adjustTimeStamps(newRegions, windows(ii).start,fs);

    if ~all(checkRegions(newRegions,fs))
        disp('bad adjusted newRegions');
        keyboard
    end

    % wait for clip to end before playing sound
    if ~params.quiet && exist('player','var')
        while strcmp(get(player,'Running'),'on'), pause(0.01); end;
    end
    
    if params.playsample
        for istat = 1:numel(newRegions)
            clip = getClip(songStruct, addPrePost(newRegions(istat),params));
            playSound(clip, params.fs, true);
            pause(0.1);
        end
    end
    
    if params.pause, pause; end
    pause(0.3);
    % adjust time of newRegions and append to stack
    regions = [regions; newRegions];
end
end

