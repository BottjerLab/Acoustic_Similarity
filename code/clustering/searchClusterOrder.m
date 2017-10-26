function [bestSeq, bestSeqFreq] = searchClusterOrder(events, links, params, varargin)

if nargin < 3; params = defaultParams; end;
params = processArgs(params,varargin{:});

% separate syllables into separate clusters,
% either clusters or number of possible clusters
distCutoff = 0.7;
clusterIdxs = cluster(links,'Cutoff', distCutoff,'criterion','Distance');
uClusters = unique(clusterIdxs);
nClusters = numel(uClusters) + 1;
transMatrix = zeros(nClusters); % let silence be the last transition possibility

% build up first-order transition matrix
nEvents = numel(events);
silenceIdx = nEvents;
silenceDuration = 2.0; % seconds, any
currSyllType = clusterIdxs(1);
seqList = currSyllType;
for ii = 1:nEvents-1 % list transition matrix between pairs of events
    nextSyllType = clusterIdxs(ii+1);
    if events(ii+1).start - events(ii).stop > silenceDuration
        transMatrix(currSyllType, silenceIdx) = ...
            transMatrix(currSyllType, silenceIdx) + 1;
        transMatrix(silenceIdx, nextSyllType) = ...
            transMatrix(silenceIdx, nextSyllType) + 1;
        seqList = [seqList silenceIdx nextSyllType];
    else
        transMatrix(currSyllType, nextSyllType) = ...
            transMatrix(currSyllType, nextSyllType) + 1;
        seqList = [seqList silenceIdx nextSyllType];
    end
    currSyllType = nextSyllType;
    
end

% convert to probabilities
nEventsPerCluster = sum(transMatrix,1);
probMatrix = transMatrix;
for ii = 1:nClusters
    probMatrix(:,ii) = probMatrix(:,ii) ./ nEventsPerCluster;
end

if params.plot
    imagesc(probMatrix);
end

%% find the most common subsequences of a certain length - 
% lookoup Algorithms on Strings, Trees and Sequences - Dan Gusfield
minimumLength = 3;
maximumArray = 1e5;
if nClusters^minimumLength > maximumArray, 
    error('TooManyPossibilities','We need to use a suffix tree...');
end;

nSeq = numel(seqList);
hashSubFreq = zeros(nClusters ^ minimumLength, 1);
% brute force: hash a subsequence to a key
for ii = 1:nSeq - minimumLength
    hashedIdx = enBase(seqList(ii:ii+minLength - 1)-1, nClusters) + 1;
    hashSubFreq(hashedIdx) = hashSubFreq(hashedIdx) + 1;
end
[sortFreq, idxFreq] = sort(hashSubFreq);

seqFreq = 10000;
ii = 1
bestSeq = zeros(10,minimumLength);
while seqFreq > 100 && ii <= 10
    bestSeq(ii,:) = unBase(idxFreq(ii) - 1) + 1
    bestSeqFreq(ii) = hashSubFreq(idxFreq(ii));
end

function ret = enBase(array, base) % 
    ret = nBase(array(1:end-1)) * base + nBase(end);
end

function arr = unBase(num, base) % 
    arr = [unBase(floor(num/base)) mod(num,base)];
end