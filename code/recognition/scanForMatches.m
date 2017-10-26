function matchEvents = scanForMatches(songStruct, event, modelSyll)
% event should be single region for now, this function scans for all
% possible matches within that event... still work in progress
% modelSyll is a Learner object (TODO: enforce w/ check)
% returns matched events within the context of the event
fs = 1/songStruct.interval;
% step 1: search in constant length

typicalLength = modelSyll.typicalLength(); % in seconds

evLength = getLength(event);
% handle case where model is longer than the event
if typicalLength > getLength(event)
    
end

tStep = 0.080; % in seconds

stepLocs = 0:tStep:(evLength - typicalLength);
nSteps = numel(stepLocs);

if nSteps == 0,
    warning('event is shorter than typical length');
    return;
end
%sweepScore = zeros(1,nSteps);
for ii = 1:nSteps
    stopStep = stepLocs(ii) + typicalLength;
    evStep(ii) = ...
        eventFromTimes(stepLocs(ii) + event.start,  ...
        stepLocs(ii) + typicalLength + event.start, fs);
end

[sweepScore,allScores] = modelSyll.fastScore(songStruct, evStep);

% step 2: find best matches according to criterion over smoothed similarity?
% graphical check - are we getting the songs?
figure(2);
subplot(211);
plotWaveform(getClip(event, songStruct), fs); xx = xlim;
subplot(212); plot(stepLocs, sweepScore, 'b-', ...
    xx, modelSyll.fastThreshold * [1 1], 'r-'); xlim(xx);

fprintf('for playback...\n');
markRegions(songStruct, event);


% step 3: do fine search
% get regions of low thresh
goodRegions = [0; (sweepScore < modelSyll.fastThreshold); 0];

% the start is probably closer to the beginning than the end of the interval, but the best bet is to select typical length + length of goodness
starts = find(diff(goodRegions) == 1); stops = find(diff(goodRegions) == -1) - 1;
%extraLen = (stops - starts) * tStep;

% do refinement
if isempty(starts) || isempty(stops),
    matchEvents = initEvents(0);
    return
end;
if any(stops > numel(stepLocs))
    stops(stops > numel(stepLocs)) = numel(stepLocs);
    disp('stops out of bounds');
end
evRoughGuesses = eventFromTimes(stepLocs(starts), stepLocs(stops) + typicalLength, fs);

[matchEvents, matchScores] = modelSyll.findBestRegions(songStruct, evRoughGuesses);

matchScores
modelSyll.threshold
isFineMatch = (matchScores < modelSyll.threshold);
fprintf('%d out of %d rough matches below fine threshold... ', ...
    sum(isFineMatch), numel(matchScores));

subplot(211);
hold on;
yy = ylim;
% place correct offset
for ii = 1:numel(matchEvents)
    if matchEvents(ii).start == matchEvents(ii).stop || matchEvents(ii).start == 0
        matchEvents(ii) = adjustTimeStamps(evRoughGuesses(ii), event.start + stepLocs(starts(ii)), fs);
    else
        matchEvents(ii) = adjustTimeStamps(matchEvents(ii), event.start + stepLocs(starts(ii)), fs);
    end
    plot(matchEvents(ii).start * [1 1],yy,'g-', matchEvents(ii).stop * [1 1],yy,'r-');    
end

fprintf('for playback...\n');
markRegions(songStruct, matchEvents);

matchEvents = matchEvents(isFineMatch);

hold off;

end