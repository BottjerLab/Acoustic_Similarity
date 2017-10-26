% construct empirical distributions

%% find all the cluster directories 

clustDirs = dir([pwd filesep 'cluster-*');
clustDirs = {clustDirs.name};
%% compile all the files that are 'clust' labeled, i.e. that have distMat
% structures and have the same clustering parameters... (need the distance parameters)

%% find the latest/largest clustering...
%% load the  all scores into continuous memory

%% find p