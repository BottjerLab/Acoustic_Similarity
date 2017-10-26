function ccmean = syllableCrosscorr(songStruct, event, syllables, params, varargin)
%% parameter handling
if nargin < 4; params = defaultParams; end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;


%% get spectrograms
[subSylls, subIdxs] = getSubEvents(event, syllables);
nSub = numel(subSylls);
for ii = 1:nSub
    fprintf('Getting spectrum for syllable %d...\n',subIdxs(ii));
    [foo,spec{ii}] = cleanClip(subSylls(ii));
end


%% get cocorrelation matrix
for ii = 1:nSub
    for jj = ii+1:nSub
        fprintf('Examining pair (%d,%d)...\n',ii,jj);
        cc{ii,jj} = similarityScan(spec{ii},spec{jj});
        cc{jj,ii} = [0.0]; % just to fill it
    end
    cc{ii,ii} = [0.0]; % just to fill it
end

%% make reductions
ccmean = cellfun(@(x) mean(x),cc);
ccmax = cellfun(@(x) max(x),cc);

% ordinalize
cclist = sort([cc{:}]);
ccOrdMean = arrayfun(@(x) find(x <= cclist,1), ccmean);
ccOrdMax = arrayfun(@(x) find(x <= cclist,1), ccmax);

% use anchor syllables
%% display
subplot(2,2,1);
imagesc(ccmean)
colorbar
title('mean')
subplot(2,2,2);
imagesc(ccmax)
colorbar
title('max')
subplot(2,2,3);
imagesc(ccOrdMean)
colorbar
title('mean ordinalize')
subplot(2,2,4);
imagesc(ccOrdMax)
colorbar
title('max ordinalize')

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