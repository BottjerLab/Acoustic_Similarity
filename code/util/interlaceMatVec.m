function matNew = interlaceMatVec(mat, vec, newVal)
% prereqs: vec2 is sparse, newVal is not zero, mat is square(?)
% vec should be a 1xN dimension, mat should be NxN
% inserts vector elements in front of/before matrix rows and columns
isNew = vec ~= 0;
nNew = sum(isNew);
origOffsets = cumsum(isNew);

newLocs = find(isNew) + (0:nNew-1);
if ~isstruct(newVal)
    matNew = zeros(size(mat,1) + nNew, size(mat,2) + nNew);
end
matNew((1:size(mat,1)) + origOffsets, (1:size(mat,2)) + origOffsets) = mat;
matNew(newLocs,:) = newVal;
matNew(:,newLocs) = newVal;

if ischar(newVal)
    matNew = char(matNew);
end
end