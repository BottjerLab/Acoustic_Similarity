function forest = trainClassifier(songStruct, events, features, marks)

% here features are the reduction of feature vectors, i.e. 
% it is expected that features is a struct array where each entry for all fields 
% is a single value
%% setup and checks
assert(numel(events) == numel(features));
if nargin < 4
    marks = markRegions(songStruct, events);
end

% convert true/false logical array to -1/1 array
if islogical(marks)
   newMarks = ones(size(marks));
   newMarks(~marks) = -1;
   marks = newMarks;
end
assert(numel(marks) == numel(features));

fields = fieldnames(features);
nF = numel(fields); % number of features

nE = numel(marks); % number of examples

%%
% some free parameters
trainFrac = 0.8; % how much of data will be training, higher means more accurate, but more overfitting
maxIters = 100;

%% compile feature matrix
featVecs = zeros(nF, nE); % columns-major
for ii = 1:nF
    featVecs(ii,:) = [features.(fields{ii})];
end

%% get training/testing splits
nTrain = floor(trainFrac * nE);
[y,ii] = sort(rand(nE,1)); % get random ordering
iTrain = ii(1:nTrain);
iTest = ii(nTrain+1:end);

% split training/testing data
trainData   = featVecs(:,iTrain);
trainLabels =      marks(iTrain);

testData   = featVecs(:,iTest);
testLabels =      marks(iTest);

%% constructing weak learner
weak_learner = tree_node_w(3); % pass the number of tree splits to the constructor

% train with Modest AdaBoost ( does this tolerate NaNs?)
[forest.Learners forest.Weights] = GentleAdaBoost(weak_learner, trainData, trainLabels, maxIters);

%% get results and make confusion matrix
fullResult = sign(Classify(forest.Learners, forest.Weights, featVecs));

fprintf('Full confusion matrix over all data (rows = true values, columns = fitted values)\n');
confMatrix = zeros(2,2);
for ii = -1:2:1
    for jj = -1:2:1
        confMatrix((ii+3)/2,(jj+3)/2) = sum(marks == ii & fullResult == jj);
    end
end
confMatrix

end