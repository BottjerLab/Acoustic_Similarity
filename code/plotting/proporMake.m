function [props1, props2, n1, n2] = proporMake(pVal, rSVal, bigBirdStageData, isCoreAndSUAFlag, isShellAndSUAFlag)
%proporMake JMA 12/11/13
%   getting together all the proportions for plotting
x = [bigBirdStageData(isCoreAndSUAFlag).(pVal)] < 0.05 & [bigBirdStageData(isCoreAndSUAFlag).(rSVal)] > 0;
props1(1) = sum(x) / sum(isCoreAndSUAFlag);
n1(1) = sum(x);
y = [bigBirdStageData(isCoreAndSUAFlag).(pVal)] < 0.05 & [bigBirdStageData(isCoreAndSUAFlag).(rSVal)] < 0;
props2(1) = sum(y) / sum(isCoreAndSUAFlag);
n2(1) = sum(y);
x = [bigBirdStageData(isShellAndSUAFlag).(pVal)] < 0.05 & [bigBirdStageData(isShellAndSUAFlag).(rSVal)] > 0;
props1(2) = sum(x) / sum(isShellAndSUAFlag); 
n1(2) = sum(x);
y = [bigBirdStageData(isShellAndSUAFlag).(pVal)] < 0.05 & [bigBirdStageData(isShellAndSUAFlag).(rSVal)] < 0;
props2(2) = sum(y) / sum(isShellAndSUAFlag);
n2(2) = sum(y);
end

