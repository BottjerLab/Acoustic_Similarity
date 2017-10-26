function spectrumData = getSpectrumStats(sample, params)
% acquire spectrogram - TODO: use multitaper
% https://github.com/dmeliza/libtfr
winSS = floor(params.windowSize * params.fs/1000);
overlapSS = floor(params.nOverlap * params.fs/1000);
if ~isfield(params,'freqBands')
    [spectrumData.spectrum, spectrumData.freqs, spectrumData.times, spectrumData.psd] = ...
        spectrogram(sample, winSS, overlapSS, params.NfreqBands, params.fs);
else
    [spectrumData.spectrum, spectrumData.freqs, spectrumData.times, spectrumData.psd] = ...
        spectrogram(sample, winSS, overlapSS, params.freqBands, params.fs);
end

% get center frequency - simply weighted average by spectrogram
spectrumData.centerFreq = getCentroidFreq(spectrumData);

% get total power
spectrumData.totalPower = getTotalPower(spectrumData);

% get wiener entropy: high = noise, low = pure tone, mid = stack/chirp ramp
spectrumData.wienerEntropy = getEntropy(spectrumData);

% get derivative spectrogram
%spectrumData.deriv = getSpectralDerivative(spectrumData, params);
[spectrumData.deriv spectrumData.mTD, spectrumData.mFD, spectrumData.FM] ...
     = getSpectralDerivative(spectrumData);

 % get pitch goodness
[spectrumData.pitchGoodness, spectrumData.harmonicPitch] = ...
    getHarmonicPitch(spectrumData, params);
end

function power = getTotalPower(spectrum)
dFreq = spectrum.freqs(2) - spectrum.freqs(1);
power = sqrt(dot(spectrum.psd,spectrum.psd)) * dFreq; % amp factor
end

function centerfreq = getCentroidFreq(spectrum)
centerfreq = spectrum.freqs' * spectrum.psd ./ sum(spectrum.psd);
end

function entropy = getEntropy(spectrum)
% Wiener entropy, defined as the log ratio of GM to AM of the power
% pure white noise means that Wiener entropy should be 0, and pure tone
% should have -inf Wiener entropy

AMpsd = mean(spectrum.psd,1);
logGMpsd = mean(log(spectrum.psd),1);
entropy = logGMpsd - log(AMpsd);
end

% get mag gradient of spectral derivative, as well as frequency Modulation
function [deriv, mTD, mFD, freqMod] = getSpectralDerivative(spectrum)
% and maximum derivatives in time and space
dt = spectrum.times(2) - spectrum.times(1);
df = spectrum.freqs(2) - spectrum.freqs(1);
freqDeriv = conv2([1 0 -1], [1], spectrum.psd, 'same') / df;
timeDeriv = conv2([1], [1 0 -1], spectrum.psd, 'same') / dt;
deriv = sqrt(timeDeriv.*timeDeriv + ...
    freqDeriv .* freqDeriv) .* sign(freqDeriv);

mTD = max(abs(timeDeriv),[],1);
mFD = max(abs(freqDeriv),[],1);
%[foo,imax] = max(abs(deriv), [], 1);
% this definition doesn't seem to be that sensitive
freqMod = atan2(max(abs(freqDeriv),[],1),max(abs(timeDeriv / 1000),[],1)) * 180 / pi;

%{
tPart = sin(freqMod);
fPart = cos(freqMod);
deriv = timeDeriv .* (ones(size(spectrum.freqs)) * tPart) + ...
    freqDeriv .* (ones(size(spectrum.freqs)) * fPart) ./ spectrum.psd;
%}
end


function [goodness, pitch] = getHarmonicPitch(spectrum, params)
dcepstrum = fft(spectrum.deriv,[],1);
dcepstrum = dcepstrum(1:ceil(end/2),:); % remove the non-symmetric part out

% pitch goodness is unscaled, calculated from deriv-cepstrum as in SAP
% 2011
% harmonic pitch is estimated as well - doesn't seem to work too well

[goodness, pitch] = max(abs(dcepstrum),[],1);
pitch = params.fs./pitch;
end

