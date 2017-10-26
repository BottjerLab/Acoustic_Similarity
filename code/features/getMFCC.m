function mfcc = getMFCC(spectrum, specparams)

% code adapted from the auditory toolbox, (c) malcolm stanley 1993
%	Filter bank parameters
lowestFrequency = 400/3;
linearFilters = 13;
linearSpacing = 200/3;
logFilters = 27;
logSpacing = 1.0711703;
fftSize = numel(spectrum.freqs);
cepstralCoefficients = 13;

windowSize = 256;		% Standard says 400, but 256 makes more sense
				% Really should be a function of the sample
				% rate (and the lowestFrequency) and the
				% frame rate.
%if (nargin < 2) 
samplingRate = specparams.fs;
windowStep = (specparams.windowSize / 1000) * samplingRate; % input is in ms, should be in samples
if (nargin < 3) frameRate = 100; end;

% Keep this around for later....
totalFilters = linearFilters + logFilters;

% Now figure the band edges.  Interesting frequencies are spaced
% by linearSpacing for a while, then go logarithmic.  First figure
% all the interesting frequencies.  Lower, center, and upper band
% edges are all consequtive interesting frequencies. 

freqs = lowestFrequency + (0:linearFilters-1)*linearSpacing;
freqs(linearFilters+1:totalFilters+2) = ...
		      freqs(linearFilters) * logSpacing.^(1:logFilters+2);

lower = freqs(1:totalFilters);
center = freqs(2:totalFilters+1);
upper = freqs(3:totalFilters+2);

% We now want to combine FFT bins so that each filter has unit
% weight, assuming a triangular weighting function.  First figure
% out the height of the triangle, then we can figure out each 
% frequencies contribution
mfccFilterWeights = zeros(totalFilters,fftSize);
triangleHeight = 2./(upper-lower);
fftFreqs = spectrum.freqs;

for chan=1:totalFilters
    mfccFilterWeights(chan,:) = ...
  (fftFreqs > lower(chan) & fftFreqs <= center(chan)).* ...
   triangleHeight(chan).*(fftFreqs-lower(chan))/(center(chan)-lower(chan)) + ...
  (fftFreqs > center(chan) & fftFreqs < upper(chan)).* ...
   triangleHeight(chan).*(upper(chan)-fftFreqs)/(upper(chan)-center(chan));
end

% try graphing
%semilogx(fftFreqs,mfccFilterWeights')
%axis([lower(1) upper(totalFilters) 0 max(max(mfccFilterWeights))])

% hamWindow = 0.54 - 0.46*cos(2*pi*(0:windowSize-1)/windowSize);

% Figure out Discrete Cosine Transform.  We want a matrix
% dct(i,j) which is totalFilters x cepstralCoefficients in size.
% The i,j component is given by 
%                cos( i * (j+0.5)/totalFilters pi )
% where we have assumed that i and j start at 0.
mfccDCTMatrix = 1/sqrt(totalFilters/2)*cos((0:(cepstralCoefficients-1))' * ...
				(2*(0:(totalFilters-1))+1) * pi/2/totalFilters);
mfccDCTMatrix(1,:) = mfccDCTMatrix(1,:) * sqrt(2)/2;

cols = numel(spectrum.times);
mfcc = zeros(cepstralCoefficients, cols);

% Ok, now let's do the processing.  For each chunk of data:
%    * Shift it into FFT order,
%    * Find the magnitude of the fft,
%    * Convert the fft data into filter bank outputs,
%    * Find the log base 10,
%    * Find the cosine transform to reduce dimensionality.
%for start=0:cols-1

%    * Window the data with a hamming window,
%    first = start*windowStep + 1;
%    last = first + windowSize-1;
%    fftData = zeros(1,fftSize);
%    fftData(1:windowSize) = preEmphasized(first:last).*hamWindow;

%    fftMag = abs(fft(fftData));
    fftMag = spectrum.psd;

% apply the triangular weights
for start = 1:cols
    earMag = log10(mfccFilterWeights * fftMag(:,start));
    mfcc(:, start) = mfccDCTMatrix * earMag;
%    mfcc(:,start+1) = mfccDCTMatrix * earMag;
end
