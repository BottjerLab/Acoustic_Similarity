for ii = 1:numel(qqq)
    [newFeats,newSylls]=getFeatures(Lb277_3_27_4_Ch1,qqq(ii),noiseProfile,...
        'plot',false,'verbose',true,'playsample',false);
    if ~isempty(newFeats)
        noisy = testIsNoise(forest, newFeats(1));
    end
end

function isNoise = testIsNoise(Forest, featureStruct)
%% compile feature matrix
fields = fieldnames(featureStruct);
nF = numel(fields); % number of features
nE = numel(marks); % number of examples
featVecs = zeros(nF, nE); % columns-major
for ii = 1:nF
    featVecs(ii,:) = [features.(fields{ii})];
end
isNoise = Classify(Forest.Learners, Forest.Weights, featVecs) > 0;
end
