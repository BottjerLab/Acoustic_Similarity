function twoSetsBar(core1, shell1, core2, shell2, label, y, stageName)
%twoSetsBar JMA 12/10/13
%   plots 2 sets of core and shell mean values with sem and ttest

m1 = [nanmean(core1) nanmean(shell1)];
SEM1 = [nanstd(core1)/sqrt(sum(~isnan(core1))) nanstd(shell1)/sqrt(sum(~isnan(shell1)))];

m2 = [nanmean(core2) nanmean(shell2)];
SEM2 = [nanstd(core2)/sqrt(sum(~isnan(core2))) nanstd(shell2)/sqrt(sum(~isnan((shell2))))];

[t1, p1] = ttest2(core1,shell1);
pVal1 = {['p = ' num2str(p1,3)]};
if ~isempty(core2) & ~isempty (shell2)
[t2,p2] = ttest2(core2,shell2);
tVals = [t1,t2];
pVal2 = {['p = ' num2str(p2,3)]};
else
    tVals = t1;
    pVal2 = 'N/A';
end

barError = [SEM1; SEM2];
barFraction = [m1; m2];
figure;
plotBarError(barFraction,barError,[],tVals,[0.7 0.7 0.7; 1 0 0]);
set(gca, 'XTickLabel', label);
ylabel(y);
legend({['core (n = ' num2str(sum(~isnan(core1))) ', ' num2str(sum(~isnan(core2))) ')'],...
    ['shell (n = ' num2str(sum(~isnan(shell1))) ', ' num2str(sum(~isnan(shell2))) ')']});
xx = get(gca,'XTick');
yy = get(gca, 'YLim');
if yy(1) == 0
    yy(1) = 0.05*yy(2);
end
text(xx(1),yy(1)+ .1*abs(yy(1)),pVal1);
text(xx(2),yy(1)+ .1*abs(yy(1)),pVal2);

title(stageName);
end

