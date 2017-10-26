function [meanCon, varCon] = conformity(clusterMat, types)

% under construction

if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end
% isSeld must be a logical

N = numel(types);
for ii = 1:N
    iType = types(ii);
    isAlsoType = types == iType;
    isNotType = ~isAlsoType;
    isAlsoType(ii) = false;
    meanCon(ii) = mean(clusterMat(ii,isAlsoType))/ mean(clusterMat(ii,isAlsoType | isNotType));
    varCon(ii)  = var (clusterMat(ii,isAlsoType))/ var (clusterMat(ii,isAlsoType | isNotType));
end