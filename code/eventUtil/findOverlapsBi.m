function [overlapIndices, overlapAreas] = findOverlapsBi(eventsA, eventsB)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% NB: eventsA and eventsB should be column vectors;

nA = numel(eventsA);
nB = numel(eventsB);

% not the most efficient but whatever
[overlapInds, overlapAreas] = findOverlaps([eventsA; eventsB]);

% pick only pairs that represent bipartite edges
overlapInds = sort(overlapInds, 1);
isBiPair = (overlapInds(:,1) <= nA & overlapInds(:,2) >  nA);           
overlapIndices = overlapInds (isBiPair,:);
overlapAreas   = overlapAreas(isBiPair  );
% unshift the offset of the second set of events
overlapIndices(:,2) = overlapIndices(:,2) - nA;
end