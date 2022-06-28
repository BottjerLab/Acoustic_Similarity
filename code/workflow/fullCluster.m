function [syllClusters] = fullCluster(songStruct,DRsylls,noiseMask,featureTable,spectra)
% Function to cluster the syllables into various groups.
%
% Requires the following inputs:
%   songStruct:
%   DRsylls:
%   noiseMask:
%   featureTable: the struct containing the features generated in
%                 prepareDRsylls.m
%   spectra: the struct containing the spectra generated in
%            prepareDRsylls.m
%
% Edited by EL 2021
%%
params = processArgs(defaultParams, 'warpingCost', 1.0);

t1 = clock;
[clusterIdxs, empMats, distMats, empDistrs] = DRcluster(songStruct,DRsylls,noiseMask,featureTable,spectra,params);
fprintf('Time for total clustering: %0.2fs\n',etime(clock, t1));

% make some figures and save work
syllClusters.clusterIdxs = clusterIdxs;
syllClusters.empMats = empMats;
syllClusters.distMats = distMats;
syllClusters.empDistrs = empDistrs;

%     save('altClust.mat','syllClusters');
%%
% takenIdxs = clusterIdxs(:,end); %
% nTypes = max(takenIdxs);
% 
% for jj = 1:nTypes
%     fig = figure(jj);
%     theseSylls = DRsylls(takenIdxs==jj);
%     fprintf('Plotting trained clusters for syllable %d (#=%d)...\n',jj, sum(takenIdxs==jj));
%     
%     mosaicDRSpec(theseSylls, songStruct, params, 'dgram.minContrast', 1e-10, ...
%         'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5, 'doFilterNoise',true,'noiseFilter', noiseMask);
%     set(fig, 'Name', sprintf('Syllable #%d, Full',jj));
%     figFileName = [clusterFile '-a' num2str(thisAge) '-c' num2str(jj) '-train.jpg'];
%     fprintf('Saving figure to %s...\n',figFileName);
%     saveCurrFigure(figFileName);
%     close(fig);
% end
end