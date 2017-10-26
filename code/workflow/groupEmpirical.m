clustDirs = dir([pwd filesep 'data' filesep 'cluster-*']);
clustDirs = {clustDirs.name};

% load all the empirical distributions
% find each file containing empirical distributions, prefixed clustDataAge
