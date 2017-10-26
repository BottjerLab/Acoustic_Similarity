function fourSetsBar(core1, shell1, core2, shell2, core3, shell3, core4, shell4, label, y, stageName)
%fourSetsBar JMA adapted from twoSetsBar 10/24/14
%   plots 4 sets of core and shell mean values with sem and ttest

m1 = [nanmean(core1) nanmean(shell1)];
SEM1 = [nanstd(core1)/sqrt(sum(~isnan(core1))) nanstd(shell1)/sqrt(sum(~isnan(shell1)))];

m2 = [nanmean(core2) nanmean(shell2)];
SEM2 = [nanstd(core2)/sqrt(sum(~isnan(core2))) nanstd(shell2)/sqrt(sum(~isnan((shell2))))];

m3 = [nanmean(core3) nanmean(shell3)];
SEM3 = [nanstd(core3)/sqrt(sum(~isnan(core3))) nanstd(shell3)/sqrt(sum(~isnan((shell3))))];

m4 = [nanmean(core4) nanmean(shell4)];
SEM4 = [nanstd(core4)/sqrt(sum(~isnan(core4))) nanstd(shell4)/sqrt(sum(~isnan((shell4))))];

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
if ~isempty(core3) & ~isempty (shell3)
[t3,p3] = ttest2(core3,shell3);
tVals = [t1,t2,t3];
pVal3 = {['p = ' num2str(p3,3)]};
end
if ~isempty(core4) & ~isempty (shell4)
[t4,p4] = ttest2(core4,shell4);
tVals = [t1,t2,t3,t4];
pVal4 = {['p = ' num2str(p4,3)]};
end

[t5, p5] = ttest2(core1,core2);
pVal5 = {['p = ' num2str(p5,3)]};
[t6, p6] = ttest2(shell1,shell2);
pVal6 = {['p = ' num2str(p6,3)]};
[t7, p7] = ttest2(core3,core4);
pVal7 = {['p = ' num2str(p7,3)]};
[t8, p8] = ttest2(shell3,shell4);
pVal8 = {['p = ' num2str(p8,3)]};


barError = [SEM1; SEM2; SEM3; SEM4];
barFraction = [m1; m2; m3; m4];
figure;
plotBarError(barFraction,barError,[],tVals,[0.7 0.7 0.7; 1 0 0]);
set(gca, 'XTickLabel', label);
ylabel(y);
legend({['core (n = ' num2str(sum(~isnan(core2))) ', ' num2str(sum(~isnan(core4))) ')'],...
    ['shell (n = ' num2str(sum(~isnan(shell2))) ', ' num2str(sum(~isnan(shell4))) ')']});
xx = get(gca,'XTick');
yy = get(gca, 'YLim');
if yy(1) == 0
    yy(1) = 0.05*yy(2);
end
text(xx(1),yy(1)+ .1*abs(yy(1)),pVal1);
text(xx(2),yy(1)+ .1*abs(yy(1)),pVal2);
text(xx(3),yy(1)+ .1*abs(yy(1)),pVal3);
text(xx(4),yy(1)+ .1*abs(yy(1)),pVal4);

text(xx(1)+ 0.5,yy(1)+ .9*abs(yy(1)),pVal5);
text(xx(1)+ 0.5,yy(1)+ .5*abs(yy(1)),pVal6);
text(xx(3)+ 0.5,yy(1)+ .9*abs(yy(1)),pVal7);
text(xx(3)+ 0.5,yy(1)+ .5*abs(yy(1)),pVal8);

title(stageName);
end

