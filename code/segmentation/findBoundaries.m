function [silentRegions, soundRegions] = findBoundaries(songStruct, region, params, varargin)

%under construction ....

% some wavelet analysis to find sound onsets 
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});
fs = 1/songStruct.interval;
params.fs = fs;
params.(params.editSpecType).fs = fs;
params.(params.editSpecType).features = ...
    [params.(params.editSpecType).features 'mfcc'];
% get waveform and spectrum data

clip = getClipAndProcess(songStruct, region, params);
spectrum = getMTSpectrumStats(clip, params.(params.editSpecType));

% % find amplitude modulation in different power bands
% % modifiable parameters
% nPowerBands = 7;
% scale = 5:8:61; 
% masterThresh = params.syllable.minPower;
% 
% nScales = numel(scale);
% powerBandBorders = logspace(log10(params.highPassFq), log10(params.lowPassFq), ...
%     nPowerBands + 1);
% thresholds = ones(nPowerBands) * masterThresh;
% 
% powerPerBand = zeros(nPowerBands, numel(spectrum.times), nScales);
% aboveT = zeros(nPowerBands, numel(spectrum.times), nScales);
% 
% if params.plot
%     figure(3);
%     subplot(211);
% end
% 
% cols = jet(nPowerBands);
% for ii = 1:nPowerBands
%     % powerBand = band pass filtering 
%     % todo: replace with MFCC power bands if necessary
%     includedBands = (powerBandBorders(ii) <= spectrum.freqs & spectrum.freqs < powerBandBorders(ii+1));
%     for jj = 1:nScales
%         % scale = amount of smoothing that the signal undergoes (similar to
%         % a low pass filter with large sidelobes)
%         % note: these low passes don't do very much
%         smoothedBandSignal = smoothSignal(sum(spectrum.psd(includedBands, :)),scale(jj));
%         powerPerBand(ii,:,jj) = smoothedBandSignal;
%         aboveT(ii,:,jj) = (smoothedBandSignal > thresholds(ii));
%         if params.plot
%             plot(spectrum.times, smoothedBandSignal - thresholds(ii), '-', 'Color', cols(ii, :));
%             hold on;
%         end
%     end
% end
% 
% % vote is sum over all scales 
% votePercent = sum(sum(aboveT,3),1) / (nScales * nPowerBands);
% 
% 
% % what fraction of votes are needed
% voteThresh = 1/(sqrt(nScales * nPowerBands));
% 
% if params.plot
%     hold off;
%     xlabel('time'); ylabel('power in band');
%     subplot(212);
%     plot(spectrum.times, votePercent, 'r-', ...
%         spectrum.times, voteThresh * ones(size(spectrum.times)), 'k-');
%     xlabel('time (s)');  ylabel('vote %'); 
% end
% keyboard
% % detect silence based on the times that bands are above above threshold
% diffAT = diff(votePercent > voteThresh);
% ducksThresh = find(diffAT == -1);
% jumpsThresh = find(diffAT ==  1);
% 
% [soundOnsets soundOffsets] = regBounds(ducksThresh, jumpsThresh, ...
%     numel(spectrum.times));
% [silOffsets silOnsets] = regBounds(jumpsThresh, ducksThresh, ...
%     numel(spectrum.times));
% 
% fSO = spectrum.times(1); % first sample offset
% 
% soundRegions = eventFromTimes(spectrum.times(soundOnsets)' - fSO, ...
%     spectrum.times(soundOffsets)', fs);
% silentRegions = eventFromTimes(spectrum.times(silOnsets)' - fSO, ...
%     spectrum.times(silOffsets)', fs);
% 
% if params.plot
%     plot(1)
%     plotAllFigures(spectrum, silentRegions, params);
%     % comparison figures
%     %{
% figure(2);
% parsedRegions = parseRegionsIntoSyllables(songStruct, region, params, ...
%     'nps.reduction',-15,'plot',true, 'dgram.minContrast',1e-9, ...
%     'syllable.comboLength',3,... % gap size in ms
%     'syllable.borderRise',1.2e-4,'pause',false);
% xx = spectrum.times([1 end]);
% coarseSounds = findPossibleSounds(spectrum, params);
% figure(3)
% plotAllFigures(spectrum, coarseSounds, params);
%     %}
% end
% 
% keyboard
% % adjust the timestamps before we leave
% soundRegions = adjustTimeStamps(soundRegions, region.start, fs);
% silentRegions = adjustTimeStamps(silentRegions, region.start, fs);
% 
% function [putStarts, putEnds] = regBounds(putStarts, putEnds, endNum)
% % boundary conditions - will fix later
% if numel(putStarts) == numel(putEnds) + 1
%     putEnds = [1 putEnds];
% elseif numel(putEnds) > numel(putStarts),
%     putStarts = [putStarts endNum];
% elseif numel(putEnds) == numel(putStarts) && ...
%         putEnds(1) > putStarts(1)
%     putEnds = [1 putEnds];
%     putStarts = [putStarts endNum];
% end
