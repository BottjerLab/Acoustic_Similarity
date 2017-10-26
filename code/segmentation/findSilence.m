function [silentRegions, soundRegions, spectrum] = findSilence(songStruct, region, params, varargin)

% yet another segmentation method for detecting silence...
% slower than stepSpectrogram but tends to be more accurate 
% not as accurate/slow as segmentSyllables, especially in finding boundaries
% (can be tuned in the future)
% in general tends to have false negatives in finding silence 
% (due to long smoothing)
% (which is to say, marks silence as sound)
% but has good false positive rate (doesn't mark sound as silence)
% i.e. is more robust to find silence (properly filtered)

if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});
fs = 1/songStruct.interval;
params.fs = fs;
params.(params.editSpecType).fs = fs;
% get waveform and spectrum data

clip = getClipAndProcess(songStruct, region, params);
spectrum = getMTSpectrumStats(clip, params.(params.editSpecType));

% find amplitude modulation in different power bands
% modifiable parameters
nPowerBands = 7;
scale = 5:8:61; 
masterThresh = params.syllable.minPower;

nScales = numel(scale);
powerBandBorders = logspace(log10(params.highPassFq), log10(params.lowPassFq), ...
    nPowerBands + 1);
thresholds = ones(nPowerBands) * masterThresh;

powerPerBand = zeros(nPowerBands, numel(spectrum.times), nScales);
aboveT = zeros(nPowerBands, numel(spectrum.times), nScales);

cols = jet(nPowerBands);
for ii = 1:nPowerBands
    % powerBand = band pass filtering 
    % todo: replace with MFCC power bands if necessary
    includedBands = (powerBandBorders(ii) <= spectrum.freqs & spectrum.freqs < powerBandBorders(ii+1));
    for jj = 1:nScales
        % scale = amount of smoothing that the signal undergoes (similar to
        % a low pass filter with large sidelobes)
        % note: these low passes don't do very much
        smoothedBandSignal = smoothSignal(sum(spectrum.psd(includedBands, :)),scale(jj));
        powerPerBand(ii,:,jj) = smoothedBandSignal;
        aboveT(ii,:,jj) = (smoothedBandSignal > thresholds(ii));
        if params.plot
            subplot(311);
            plot(spectrum.times, smoothedBandSignal - thresholds(ii), '-', 'Color', cols(ii, :));
            hold on;
        end
    end
end

% vote is sum over all scales 
votePercent = sum(sum(aboveT,3),1) / (nScales * nPowerBands);


% what fraction of votes are needed
voteThresh = 1/(sqrt(nScales * nPowerBands));

if params.plot %&& params.verbose
     hold off;
     xlabel('time'); ylabel('power in band');
     xlim([0 max(spectrum.times)]);
     subplot(312);
     plot(spectrum.times, votePercent, 'r-', ...
         spectrum.times, voteThresh * ones(size(spectrum.times)), 'k-');
     xlabel('time (s)');  ylabel('vote %'); 
     xlim([0 max(spectrum.times)]);
     subplot(313);
     plot((1:numel(spectrum.waveform)) / fs, spectrum.waveform, 'b-');
     xlim([0 numel(spectrum.waveform) / fs])
end

% detect silence based on the times that bands are above above threshold
diffAT = diff(votePercent > voteThresh);
ducksThresh = find(diffAT == -1);
jumpsThresh = find(diffAT ==  1);

[soundOnsets soundOffsets] = regBounds(ducksThresh, jumpsThresh, ...
    numel(spectrum.times));
[silOffsets silOnsets] = regBounds(jumpsThresh, ducksThresh, ...
    numel(spectrum.times));

fSO = spectrum.times(1); % first sample offset

soundRegions = eventFromTimes(spectrum.times(soundOnsets)' - fSO, ...
    spectrum.times(soundOffsets)', fs);
silentRegions = eventFromTimes(spectrum.times(silOnsets)' - fSO, ...
    spectrum.times(silOffsets)', fs);

% adjust the timestamps before we leave
soundRegions = adjustTimeStamps(soundRegions, region.start, fs);
silentRegions = adjustTimeStamps(silentRegions, region.start, fs);

function [putStarts, putEnds] = regBounds(putStarts, putEnds, endNum)
if numel(putStarts) == 0 && numel(putEnds) == 0, return; end
% boundary conditions - will fix later
if numel(putStarts) == numel(putEnds) + 1
    putEnds = [1 putEnds];
elseif numel(putEnds) > numel(putStarts),
    putStarts = [putStarts endNum];
elseif numel(putEnds) == numel(putStarts) && ...
        putEnds(1) > putStarts(1)
    putEnds = [1 putEnds];
    putStarts = [putStarts endNum];
end
