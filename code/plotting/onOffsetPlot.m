function onOffsetPlot(ERAOnsetWindow, onsetBinSize, RSOn, RSOff, isCoreAndSUAFlag, isShellAndSUAFlag, bigBirdStageData, type, stageName, name, sub1, sub2)
%onOffset Plot
%   plots mean onset and offset responses
ERAbins = (ERAOnsetWindow(1):onsetBinSize:ERAOnsetWindow(2))' * 1000;
onPos = RSOn > 0;
offPos = RSOff > 0;
onPosC = isCoreAndSUAFlag & onPos & sub1;
onPosS = isShellAndSUAFlag & onPos & sub1;
onNegC = isCoreAndSUAFlag & ~onPos & sub1;
onNegS = isShellAndSUAFlag & ~onPos & sub1;
offPosC = isCoreAndSUAFlag & offPos & sub2;
offPosS = isShellAndSUAFlag & offPos & sub2;
offNegC = isCoreAndSUAFlag & ~offPos & sub2;
offNegS = isShellAndSUAFlag & ~offPos & sub2;

if strcmp(type,'all')
    figure;
    subplot(211);
    plotFRERA(bigBirdStageData(onPosC),bigBirdStageData(onPosS),'ERAOnsethist',...
        ERAbins);
    xlabel(['time relative to onset (ms), positive pre-onset responses, ' name]);
    title([stageName, ' core n = ', num2str(sum(onPosC)), ' and ',num2str(sum(offPosC)),' shell n = ', num2str(sum(onPosS)),' and ', num2str(sum(offPosS))]);
    
    subplot(212);
    plotFRERA(bigBirdStageData(offPosC),bigBirdStageData(offPosS),'ERAOffsethist',...
        ERAbins);
    xlabel(['time relative to offset (ms), positive pre-offset responses ' name]);
    legend('core','shell');
    
    figure;
    subplot(211);
    plotFRERA(bigBirdStageData(onNegC),bigBirdStageData(onNegS),'ERAOnsethist',...
        ERAbins);
    xlabel(['time relative to onset (ms), negative pre-onset responses ' name]);
    title([stageName, ' core n = ', num2str(sum(onNegC)), ' and ',num2str(sum(offNegC)),' shell n = ', num2str(sum(onNegS)),' and ', num2str(sum(offNegS))]);
    
    subplot(212);
    plotFRERA(bigBirdStageData(offNegC),bigBirdStageData(offNegS),'ERAOffsethist',...
        ERAbins);
    xlabel(['time relative to offset (ms), negative pre-offset responses ' name]);
    legend('core','shell');
else
    figure;
    subplot(211);
    plotFRERA(bigBirdStageData(onPosC),bigBirdStageData(onPosS),'firstERAOnsethist',...
        ERAbins);
    xlabel(['time relative to onset (ms), positive pre-onset responses for first syllables ' name]);
    title([stageName, ' core n = ', num2str(sum(onPosC)), ' and ',num2str(sum(offPosC)),' shell n = ', num2str(sum(onPosS)),' and ', num2str(sum(offPosS))]);
    
    subplot(212);
    plotFRERA(bigBirdStageData(offPosC),bigBirdStageData(offPosS),'lastERAOffsethist',...
        ERAbins);
    xlabel(['time relative to offset (ms), positive pre-offset responses to last syllables ' name]);
    legend('core','shell');
    
    figure;
    subplot(211);
    plotFRERA(bigBirdStageData(onNegC),bigBirdStageData(onNegS),'firstERAOnsethist',...
        ERAbins);
    xlabel(['time relative to onset (ms), negative pre-onset responses for first syllables ' name]);
    title([stageName, ' core n = ', num2str(sum(onNegC)), ' and ',num2str(sum(offNegC)),' shell n = ', num2str(sum(onNegS)),' and ', num2str(sum(offNegS))]);
    
    subplot(212);
    plotFRERA(bigBirdStageData(offNegC),bigBirdStageData(offNegS),'lastERAOffsethist',...
        ERAbins);
    xlabel(['time relative to offset (ms), negative pre-offset responses to last syllables ' name]);
    legend('core','shell');
end
end

