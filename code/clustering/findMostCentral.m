function [centralOrder, subMats, avgDist] = findMostCentral(clusterMat, clusterIds)

if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end
% clusterIds must be numeric

nTypes = nanmax(clusterIds);
centralOrder = cell(1,nTypes);
subMats = cell(1,nTypes);
avgDist = cell(1,nTypes);
for ii = 1:nTypes
	thisType = find(clusterIds == ii);
	if ~any(thisType), continue; end
    
    subMats{ii} = clusterMat(thisType, thisType);
    
	% find the average distance of every element to its centroid
    avgDist{ii} = mean(subMats{ii},1);
	[~,centralOrder{ii}] = sort(avgDist{ii});
    centralOrder{ii} = thisType(centralOrder{ii});
    subMats{ii} = squareform(subMats{ii});
end