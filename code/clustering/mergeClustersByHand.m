function redoneLabels = mergeClustersByHand(syllSet, dists, currLabels)
% currLabels should have some clusters and some NaNs
plotParams = processArgs(defaultParams,...
                            'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
                            'preroll', 3, 'postroll', 3);
        
% flatten the labels to 1-N, and NaNs
nanTag = any(isnan(currLabels));
currLabels(isnan(currLabels)) = Inf;
[uLabels,~,currLabels] = unique(currLabels);
if nanTag    
    currLabels(currLabels == numel(uLabels)) = NaN;
end

nLabels = nanmax(currLabels); 
mostCentral = findMostCentral(dists, currLabels); % cell array

tStamp = datestr(clock,'mm_dd_HH_MM_SS');
% create temporary figures
fExemplar = figure;
tmpDir = [pwd filesep 'tmp-' tStamp filesep];
mkdir(tmpDir);
for ii = 1:nLabels
    nEx = min(numel(mostCentral{ii}),5);
    bestFive = syllSet(mostCentral{ii}(1:nEx));
    mosaicDRSpec(bestFive, plotParams, 'doFilterNoise', true);
    
    filName = [tmpDir 'clusterEx-' num2str(ii,'%02d') '.jpg'];
    saveCurrFigure(filName);
    clf
end    
close(fExemplar)
% use the file system to hold the figures
system(['explorer ' tmpDir]);

%%
defaultAns = cellstr(strcat(strcat('[',num2str((1:nLabels)')),'],')); 
defaultAns = ['{' strcat(defaultAns{:})]; 
defaultAns(end) = '}';

groupingStr = nm_inputdlg('Give in cell of vectors form the groupings', ...
    'Merge step',1,{defaultAns}); groupingStr = groupingStr{1};
groupings = eval(groupingStr);

fprintf('Reassigning and regrouping...')
redoneLabels = NaN(size(currLabels));
for ii = 1:numel(groupings)
    for jj = 1:numel(groupings{ii})
        redoneLabels(currLabels == groupings{ii}(jj)) = ii;
    end
end

% delete the temporary figure files
delete([tmpDir '*.jpg']); 
deleteStat = rmdir(tmpDir); % assignment to catch any error

end

