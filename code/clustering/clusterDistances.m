function [intracluster, intercluster] = clusterDistances(clusterMat, clusterIds)

% the output arguments intracluster and intercluster distances 
% are of the same size 

if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end

[clTypes, ~, rIdx] = unique(clusterIds);
nTypes = numel(clTypes);

intracluster = zeros(size(clusterIds));
intercluster = zeros(size(clusterIds));
for ii = 1:nTypes
	thisType = (rIdx == ii);
	
    intraMat = clusterMat(thisType,  thisType);
    interMat = clusterMat(thisType, ~thisType);
    
	% find the average distance of every element to its centroid
	intracluster(thisType) = mean(intraMat,1);
    outercluster = mean(interMat, 2);
    if all(size(intracluster(thisType)) == size(outercluster))
        intercluster(thisType) = intracluster(thisType) ./ outercluster;
    else
        intercluster(thisType) = intracluster(thisType) ./ outercluster';
    end
end