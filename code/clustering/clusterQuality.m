function [DB, d_i] = clusterQuality(clusterMat, clusterIds)
% DB returns the davies-bouldin index for the entire clustering scheme
% d_i returns the individual cluster qualities
if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end
% clusterIds must be numeric, not all NaN

nTypes = max(clusterIds);

centerIdx = zeros(1,nTypes);
sigma = zeros(1,nTypes);
for ii = 1:nTypes
	thisType = find(clusterIds == ii);
	subMat = clusterMat(thisType, thisType);

	% find the average distance of every element to its centroid
	%[sigma(ii), minInType] = min(mean(subMat,1));
    
    % find the RMS distance of every element to its centroid
    % NB: the centroid is defined as the element for which the rms distance 
    % from every other element is minimal
	[sigma(ii), minInType] = min(mean(subMat.^2,1));
    sigma(ii) = sqrt(sigma(ii));
    
	% find the most central element, which will be surrogate for the 'centroid'
	centerIdx(ii) = thisType(minInType);
end

% get the centroid distances
centD = clusterMat(centerIdx, centerIdx);
	
% find the d_is
d_i = zeros(1,nTypes);
for ii = 1:nTypes
	notI = 1:nTypes; notI(ii) = [];
	d_i(ii) = mean((sigma(ii) + sigma(notI)) ./ centD(ii, notI));
end
DB = mean(d_i);
end
