function [empScore, empDistr] = makeSelfEmpirical(arr, nPts)
[sArr,rord] = sort(arr);
empScore(rord) = [1:numel(arr)] / numel(arr);

xx = linspace(0,1,nPts);
if nPts == numel(arr)
    yy = sArr;
else
    yy = interp1(linspace(0,1,numel(arr)), sArr, xx);
end
% prepare for interp1 by removing redundant entries
redun = [diff(yy)==0 false];
if any(redun)
    xx = xx(~redun); % might be better to take a mean of the p-values instead of the max (as this implies)
    yy = yy(~redun);
end
empDistr = zeros(2,numel(xx));
empDistr(1,:) = xx;
empDistr(2,:) = yy;