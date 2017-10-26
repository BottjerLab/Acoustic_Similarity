function res = smoothSignal(input, filtSize, varargin)

% function res = smoothSignal(inputSize, filtSize)
% runs a gaussian smoothing window on signal INPUTSIZE, similar to low passing
% gaussian window has total width FILTSIZE, NOT SD of with FILTSIZE
% to have FILTSIZE represent standard deviations, 
% run smoothSignal(input, sdValue, 'SD');
%
% for the windowed / non-'SD' option, filtSize must be integer
% for the 'SD' option, filtSize can be non-integer.

if nargin < 2, filtSize = 7; end;

smoothType = 'nSamples';
if nargin >= 3,
    smoothType = varargin{1};
end

if filtSize == 1, res = input; return; end

binFilt = zeros(1,filtSize);
if filtSize < 20 && strcmpi(smoothType, 'nSamples')
    for i = 0:filtSize-1;
        binFilt(i+1) = nchoosek(filtSize-1,i);
    end
elseif filtSize >= 20 || strcmpi(smoothType, 'SD')
    windowSize = 128; alpha = 2 * windowSize / filtSize;
    binFilt = gausswin(windowSize, alpha); % sigma = windowSize / 2 * alpha
end
binFilt = binFilt / sum(binFilt);

% just approximate with gaussian
if versionNumber > 7.8 || versionNumber > 7.1 && versionNumber < 7.2
    padInput = input([ones(1,filtSize) 1:numel(input) numel(input) * ones(1,filtSize)]);
    
    padres = conv(padInput, binFilt, 'same');
    res = padres((filtSize+1):(end-filtSize));
else
    res = conv(input, binFilt);
    lb = length(binFilt);
    res = res((0:length(input)-1) + floor(lb/2)); % return middle
end
end