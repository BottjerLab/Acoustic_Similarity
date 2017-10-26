function [diameter, meanSep, stdSep] = clusterDiameter(clusterMat, labels, types)

% under construction
% labels must be the size of clusterMat when squared
if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end
% isSeld must be a logical

N = numel(types);
diameter = Inf(1,N);
meanSep = zeros(1,N);
stdSep = zeros(1,N);

for ii = 1:N
    thisType = types(ii);
    isType = (labels == thisType);
    
    if sum(isType) <= 1, continue; end;
        
    innerDists = squareform(clusterMat(isType, isType));
    diameter(ii) =  max(innerDists);
    meanSep(ii)  = mean(innerDists);
    stdSep(ii)   =  std(innerDists);
end