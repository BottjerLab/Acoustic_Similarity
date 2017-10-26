function forest = trainOnSubset(songStruct, events, features, trainFrac)
if nargin < 4, trainFrac = 0.25; end;
nEv = numel(events);
nSel = ceil(min(trainFrac,1) * nEv);
trainGroup = randperm(nEv,nSel);
forest = trainClassifier(songStruct, events(trainGroup), features(trainGroup));
