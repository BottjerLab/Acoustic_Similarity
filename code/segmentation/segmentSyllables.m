function regions = segmentSyllables(spectrum, params, varargin)
%SEGMENTSYLLABLES  segments spectrum objects into regions of outlying power 
% function regions = segmentSyllables(spectrum, params, varargin)

% similar to finding regions exceeding minimal power, but this function separates within the
% vocalization level, trying to identify variations in tone

% currently the function works simply by finding minima/maxima that are more extreme than
% their neighbors (interest point detection)
%
% this function all filters according to a power threshold, and a length
% threshold, and finds boundaries based on minimum increases/decreases in
% power (think of defining the borders of a mountain by where its slope
% decreases)

% in the future we can find interest points in multiple feature vectors

% assumes the spectrum is found with the .fine parameters and is
% highpassed to remove low freq noise


%% process
if nargin < 2, params = defaultParams; end;
params = processArgs(params, varargin{:});
fs = params.fs;
regions = initEvents;
%% smooth the power with simple binomial (approx gaussian) window
smoothedP = smoothSignal(spectrum.totalPower,params.syllable.smoothingWindow);

%% find minimums and quality of minimums

% look for minima/maxima - minima will become starting/ending points for
% syllables
dt = diff(spectrum.times(1:2));
diffP = [0 diff(smoothedP)]/dt;
zeroCrossings = [-diff(sign(diffP))/2 0]; % minima are negative, maxima are positive
% get the distance between each minima and its next maxima
minimaFirst = (diffP(2) < 0);  % if the first index will be a minima,i.e. if we begin by decreasing
minima = find(zeroCrossings == -1);
maxima = find(zeroCrossings ==  1);
if ~minimaFirst
    minima = [1 minima];
end
if numel(maxima) < numel(minima)
    maxima(end+1:numel(minima)) = numel(spectrum.times);
elseif numel(maxima) > numel(minima)
    maxima = maxima(1:numel(minima));
end

%% find any nearly flat minima/maxima and remove them to properly space intervals
% near flat minima/maxima are within

[minMaxSeq, sortIdx] = sort([minima maxima]);
seqIdxIsMinima = (sortIdx <= numel(minima));
dFactor = exp(abs(diff(log(smoothedP(minMaxSeq)))));
isFlat = [(dFactor < params.syllable.flatFactor) false];
doRemove = false(1,numel(minMaxSeq));
if any(isFlat)
    % remove all the extrema from sequences with an even number of flat
    % extrema
    % and leave one extrema from sequences with an odd number of flat
    % extrema
    
    pastTrue = isFlat(1);
    for ii = 2:numel(isFlat)
        if isFlat(ii), pastTrue = pastTrue + 1;
        else
            flatRange = (ii-pastTrue):ii;
            if mod(pastTrue,2) == 0;
                % is this a maxima or minima?
                isMaxima = sum(seqIdxIsMinima(flatRange)) < ...
                    sum(~seqIdxIsMinima(flatRange));
                % pick the most extreme extrema to represent the group
                
                if isMaxima
                    [foo,bestExtremaIdx] = max(smoothedP(minMaxSeq(flatRange)));
                else
                    [foo,bestExtremaIdx] = min(smoothedP(minMaxSeq(flatRange)));
                end
                doRemove(flatRange) = true;
                doRemove(flatRange(bestExtremaIdx)) = false;
            else
                doRemove(flatRange) = true;
            end
            pastTrue = 0;
        end
    end
    % make sure that an equal amount of maxima and minima are taken out
    if sum(seqIdxIsMinima(doRemove)) ~= sum(~seqIdxIsMinima(doRemove))
        fprintf('Error in reducing flat minima...\n');
        error('error','Error in reducing flat minima');
    end
    
    minMaxSeq(doRemove) = [];
    seqIdxIsMinima(doRemove) = [];
    minima = minMaxSeq(seqIdxIsMinima);
    maxima = minMaxSeq(~seqIdxIsMinima);
end


%% what happens if there are no minima?
if numel(maxima) == 0 || numel(minima) == 0
    risesThresh = []; fallsThresh = []; regions = []; 
    if params.plot,
        plotAtExit; 
    end
    return;
end

%% get length and quality
durationToMaxima = spectrum.times(maxima) - spectrum.times(minima);
QminimaFwd = smoothedP(maxima) ./ smoothedP(minima);
QminimaBkwd = smoothedP(maxima([1 1:end-1])) ./ smoothedP(minima);
Qminima = sqrt(QminimaFwd .* QminimaBkwd);

tooShortRise = (durationToMaxima < params.syllable.minRiseLength/1000);
tooQuietRise = (Qminima < params.syllable.minRiseAmp);
isBad = tooShortRise | tooQuietRise;
risesThresh = minima(~isBad);

%% set the trailing edge for each candidate syllable at the successive
% minimums
if ~isempty(risesThresh)
    fallsThresh = [risesThresh(2:end) numel(spectrum.times)];
else
    fallsThresh = zeros(1,0);
end

if numel(risesThresh) ~= numel(fallsThresh)
    fprintf('Help!');
end

%% if the syllable is too quiet, exclude it, (according to minPower)
tooSoft = false(1,numel(risesThresh));
for ii = 1:numel(risesThresh)
    tooSoft(ii) = max(smoothedP(risesThresh(ii):fallsThresh(ii))) < params.syllable.minPower;
end
risesThresh(tooSoft) = []; fallsThresh(tooSoft) = [];

risesDraft = risesThresh; fallsDraft = fallsThresh;
%% debugging - assess minima finding 
if params.plot
    figure(2)
    semilogy(1:numel(smoothedP), smoothedP, 'b-', ...
        risesThresh, smoothedP(risesThresh) * 3, 'g*',fallsThresh,smoothedP(fallsThresh) * 3.3,'r*');
    xlim([1 numel(smoothedP)])
    figure(1)
end
%% contract syllables to major rises/dropoffs, as defined by syllable.borderRise
%unsmoothedDiffP = [0 diff(spectrum.totalPower)] / dt;
smoothedDiffP = [0 diff(smoothedP)] / dt;
bestRising = (smoothedDiffP > params.syllable.borderRise);
bestFalling = (smoothedDiffP < -params.syllable.borderRise);
range = 1:numel(spectrum.times);

if numel(risesThresh) ~= numel(fallsThresh)
    error('unequal number of interest points at start and stop');
end
for ii = 1:numel(risesThresh)
    
    newRise = find(bestRising(risesThresh(ii):end),1) + risesThresh(ii);
    newFall = find(bestFalling & range <= fallsThresh(ii),1,'last');
    
    % contract the end first, but enforce order    
    if ~isempty(newFall) && newFall > risesThresh(ii), fallsThresh(ii) = newFall; end;
    if ~isempty(newRise) && newRise < fallsThresh(ii), risesThresh(ii) = newRise; end;
end

risesEdit = risesThresh; fallsEdit = fallsThresh;
%% check to make sure the number of start and end markers are consistent
if numel(risesThresh) ~= numel(fallsThresh)
    error('segmentSyllables:inconsistentRegions',...
        'Number of end markers not equal to number of start markers');
end

% if we have no regions, just quit here
if isempty(risesThresh),
    regions = initEvents;
    return;
end

% merge if there are adjacent regions (identical samples)
idxsStuckAtEnd = find(fallsThresh(1:end-1) == risesThresh(2:end));
if ~isempty(idxsStuckAtEnd)
    fallsThresh(idxsStuckAtEnd) = []; 
    risesThresh(idxsStuckAtEnd + 1) = []; 
end

%% more debugging
if params.plot
    figure(2)
    hold on;
    semilogy(risesThresh , smoothedP(risesThresh) * 1.8, 'g^',fallsThresh,smoothedP(fallsThresh) * 1.2,'r^');
    hold off;
    figure(1)
end
%% removing syllables that might be too short
%{
finalLengths = spectrum.times(fallsThresh) - spectrum.times(risesThresh);
tooShort = (finalLengths < params.syllable.minLength/1000);
if params.verbose && any(tooShort)
    fprintf('Removing %d regions that are too short ...\n', sum(tooShort));
end
fallsThresh(tooShort) = []; risesThresh(tooShort) = [];
%}
% initialize region structure
regions = eventFromTimes(spectrum.times(risesThresh), spectrum.times(fallsThresh), fs);

if params.plot 
    plotAtExit;
end

    function plotAtExit
        [hax, optGs] = plotAllFigures(spectrum, regions, params);
        
        % plot threshold for syllable detection
        axes(hax(strcmp('totalPower', optGs)));
        hold on;
        semilogy(xlim, params.syllable.minPower * [1 1],'b--');
        semilogx(spectrum.times(risesThresh), spectrum.totalPower(fallsThresh), 'kv','MarkerSize',5);
        semilogx(spectrum.times(risesThresh), spectrum.totalPower(fallsThresh), 'g^','MarkerSize',5);
        hold off;
        axes(hax(1));
    end
end