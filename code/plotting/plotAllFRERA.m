function plotAllFRERA(ERAOnsetWindow, onsetBinSize, RSOn, RSOff, isCoreAndSUAFlag, isShellAndSUAFlag, bigBirdStageData, stageName, name, sub1, sub2)
% plotAllFRERA JMA 12/11/13
%plots all the onset/offset firing rates

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

figure;
subplot(211);
on = [bigBirdStageData(onPosC).ERAOnsethist];
for x = 1: size(on,2)
    plot(ERAbins(1:end-3), on(:,x), '-', 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlim([-150 150]);
xlabel(['time relative to onset (ms), positive pre-onset responses, ' name]);
title([stageName, ' core n = ', num2str(sum(onPosC)),' shell n = ', num2str(sum(onPosS))]);
subplot(212);
on = [bigBirdStageData(onPosS).ERAOnsethist];
for y = 1:size(on,2)
    plot(ERAbins(1:end-3), on(:,y), '-', 'Color', [1 0 0]);
    hold on;
end
xlabel(['time relative to onset (ms), positive pre-onset responses, ' name]);
xlim([-150 150]);
hold off;

figure;
subplot(211);
off = [bigBirdStageData(offPosC).ERAOffsethist];
for x = 1: size(off,2)
    plot(ERAbins(1:end-3),off(:,x), '-', 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlim([-150 150]);
xlabel(['time relative to onset (ms), positive pre-offset responses, ' name]);
title([stageName, ' core n = ', num2str(sum(offPosC)),' shell n = ', num2str(sum(offPosS))]);
subplot(212);
off = [bigBirdStageData(offPosS).ERAOffsethist];
for y = 1:size(off,2)
    plot(ERAbins(1:end-3), off(:,y), '-', 'Color', [1 0 0]);
    hold on;
end
xlabel(['time relative to onset (ms), positive pre-offset responses, ' name]);
xlim([-150 150]);
hold off;

figure;
subplot(211);
on = [bigBirdStageData(onNegC).ERAOnsethist];
for x = 1: size(on,2)
    plot(ERAbins(1:end-3), on(:,x), '-', 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlim([-150 150]);
xlabel(['time relative to onset (ms), negative pre-onset responses, ' name]);
title([stageName, ' core n = ', num2str(sum(onNegC)),' shell n = ', num2str(sum(onNegS))]);
subplot(212);
on = [bigBirdStageData(onNegS).ERAOnsethist];
for y = 1:size(on,2)
    plot(ERAbins(1:end-3), on(:,y), '-', 'Color', [1 0 0]);
    hold on;
end
xlabel(['time relative to onset (ms), negative pre-onset responses, ' name]);
xlim([-150 150]);
hold off;

figure;
subplot(211);
off = [bigBirdStageData(offNegC).ERAOffsethist];
for x = 1: size(off,2)
    plot(ERAbins(1:end-3),off(:,x), '-', 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlim([-150 150]);
xlabel(['time relative to onset (ms), negative pre-offset responses, ' name]);
title([stageName, ' core n = ', num2str(sum(offNegC)),' shell n = ', num2str(sum(offNegS))]);
subplot(212);
off = [bigBirdStageData(offNegS).ERAOffsethist];
for y = 1:size(off,2)
    plot(ERAbins(1:end-3), off(:,y), '-', 'Color', [1 0 0]);
    hold on;
end
xlabel(['time relative to onset (ms), negative pre-offset responses, ' name]);
xlim([-150 150]);
hold off;
end


