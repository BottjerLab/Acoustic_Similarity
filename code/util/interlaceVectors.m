function vecNew = interlaceVectors(vec1, vec2, newVal)
% prereqs, vec2 is sparse, newVal is not zero

isNew = vec2 ~= 0;
nNew = sum(isNew);
origOffsets = cumsum(isNew);

newLocs = find(isNew) + (0:nNew-1);
if ~isstruct(newVal)
    vecNew = zeros(1, numel(vec1) + nNew);
end
vecNew((1:numel(vec1)) + origOffsets) = vec1;
vecNew(newLocs) = newVal;
if ischar(newVal)
    vecNew = char(vecNew);
end
end