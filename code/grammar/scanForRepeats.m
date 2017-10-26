function sim = scanForRepeats(songStruct, events, varargin)
%% parameter handling
if nargin < 5; params = defaultParams; end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;

%% preprocessing
% filter out noise first
[firstClip,firstSpec] = cleanClip(events(1));

for ii = 2:numel(events)
    fprintf('Checking similarity between regions %d and %d...\n', ii-1, ii);
    [nextClip,nextSpec] = cleanClip(events(ii));
    sim{ii-1} = similarityScan(firstSpec, nextSpec);
    firstSpec = nextSpec;
end

    function [clip,spec] = cleanClip(region)
        clip = getClip(addPrePost(region,params), songStruct);       
        
        if ~params.quiet, playSound(clip,fs); end
        
        % high pass the result
        params.fs=fs;
        clip = highPassSample(clip,params);
        
        %% spectral analysis
        params.fine.fs = fs;
        spec = getMTSpectrumStats(clip, params.fine);
    end
end