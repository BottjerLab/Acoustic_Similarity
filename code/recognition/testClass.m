function [isOfType, evidenceValue] = testClass(Forest, featureStruct)
%% compile feature matrix
fields = fieldnames(featureStruct);
nF = numel(fields); % number of features
nE = numel(featureStruct); % number of examples
featVecs = zeros(nF, nE); % columns-major
for ii = 1:nF
    featVecs(ii,:) = [featureStruct.(fields{ii})];
end
evidenceValue = Classify(Forest.Learners, Forest.Weights, featVecs);
isOfType = evidenceValue > 0;
end
