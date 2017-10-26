function regions = findPossibleSounds(spectrum, params, varargin)
% FINDPOSSIBLESOUNDS segment sounds roughly
%   findPossibleSounds(spectrum) segments sounds based on raw amplitude  
% and fast amplitude changes in single band
% 
% Currently, three criteria are used:
% 1a) (power must be higher than a threshold OR
% 1b) change in power must be higher than a threshold ) AND
% 2) region is not too small and 
% 3) centroid is frequency is not too low (almost always cage noise)
% Furthermore, sounds are extended to where sounds reach a minimum taper
% (i.e. they must settle to a minimum rise/fall to be cut off.)

if nargin < 2
    params = defaultParams;
end
if nargin > 2
    params = processArgs(params, varargin{:});
end

% regions begins life as a struct array
regions = initEvents;

% smooth the power with simple binomial (approx gaussian) window
smoothedP = smoothSignal(spectrum.totalPower,5);

dt = spectrum.times(2) - spectrum.times(1);
dPdT = [diff(smoothedP) 0] / dt;

% use marking based on simple derivative threshold
exceedsDiffThresh = (dPdT > params.riseThresh);
onsets = find([diff(exceedsDiffThresh) 0] == 1);
offsets = zeros(1,numel(onsets));

% for each onset, find the place where the power returns to its original
for ii = 1:numel(onsets)
    len = find(smoothedP(onsets(ii)+1:end) < smoothedP(onsets(ii)),1);
    if ~isempty(len)
        offsets(ii) = onsets(ii) + len;
    else
        offsets(ii) = numel(spectrum.times);
    end
end

% mark based on simple power threshold
exceedsPowerThresh = (smoothedP > params.powerThresh);
crossesThresh = [diff(exceedsPowerThresh) 0];
risesThresh = find(crossesThresh == 1);
fallsThresh = find(crossesThresh == -1);

% correct power thresholds
if numel(risesThresh) + numel(fallsThresh) > 0
    % boundary conditions - will fix later
    if numel(fallsThresh) == numel(risesThresh) + 1
        risesThresh = [1 risesThresh];
    elseif numel(risesThresh) > numel(fallsThresh),
        fallsThresh = [fallsThresh numel(spectrum.times)];
    elseif numel(risesThresh) == numel(fallsThresh) && ...
            risesThresh(1) > fallsThresh(1)
        risesThresh = [1 risesThresh];
        fallsThresh = [fallsThresh numel(spectrum.times)];
    end
end

% now merge the two derivative and power requirements
statChange = histc([risesThresh onsets], 1:numel(spectrum.times)) - ...
    histc([fallsThresh offsets], 1:numel(spectrum.times));
embeddedState = cumsum(statChange);
embeddedState(embeddedState>1) = 1;
if any(embeddedState) < 0 || embeddedState(end) ~= 0, error('findPossibleSounds:BadIntervals','intervals are not nested properly'); end;

onsets = find(diff([0 embeddedState]) == 1);
offsets = 1 + find(diff(embeddedState) == -1);


%%%%%%%%%%%%%%%%%% rejection code %%%%%%%%%%%%%%%%%%%%%%%
% reject any samples that are too short
minSamples = ceil(params.minLength / (params.rough.windowSize - params.rough.nOverlap));
if(numel(offsets) ~= numel(onsets)), keyboard; end;
isTooShort = (offsets - onsets < minSamples);

% reject any samples whose centroid frequencies are too low
isTooLow = false(1,numel(onsets));
for jj = 1:numel(onsets)
    isTooLow(jj) = mean(spectrum.centerFreq(onsets(jj):offsets(jj)) < ...
        params.minCenterFreq);
end
fprintf('%d region%s too short, %d region%s too low\n', sum(isTooShort), 's'*(sum(isTooShort)~=1),...
    sum(isTooLow),'s'*(sum(isTooLow)~=1));

% cut out the rejected segments
onsets(isTooShort | isTooLow)=[];
offsets(isTooShort | isTooLow)=[];

% extend all syllables to where the power is increasing
dPower = diff(smoothedP) / (params.rough.windowSize - params.rough.nOverlap);
for istat = 1:numel(onsets)
    newRise = find([ -1 dPower(1:onsets(istat))] < params.minRise, ...
        1, 'last');
    newFall = find([dPower(offsets(istat):end) 1 ] > -params.minRise, ...
        1, 'first') + offsets(istat) - 1;
    onsets(istat) = newRise;
    offsets(istat) = newFall;
end

regions = eventFromTimes(spectrum.times(onsets), spectrum.times(offsets), params.fs);
if params.verbose
    fprintf('%d regions added\n', numel(regions));
end

end
