function reducedStatTable(statNames, testStatistic)
% use fieldnames to get statNames
[statRoot, statReduction] = strtok(statNames,'_');
statReduction = strrep(statReduction,'_','');

isExcludedStat = cellfun('isempty', statReduction);

testStatistic = testStatistic(~isExcludedStat);
statRoot = statRoot(~isExcludedStat);
statReduction = statReduction(~isExcludedStat);

[uRoot, ~, iURoot] = unique(statRoot);  % statRoot = uRoot(iURoot)
[uRed , ~, iURed ] = unique(statReduction);

nURoot = numel(uRoot); nURed = numel(uRed);

% put into the grid for plotting - indexing sleight of hand
sigGrid = zeros(nURoot, nURed);
sigGrid(sub2ind([nURoot, nURed], iURoot, iURed)) = testStatistic;

% do the plotting
imagetext(sigGrid,[],[],get(gca,'Position'),false)
%colormap('summer'); % c o l o r s
xlabel('Reduction');
ylabel('Feature');
set(gca,'XTickLabel', uRed, 'YTickLabel', uRoot);
