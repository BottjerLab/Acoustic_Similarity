function corrMatrix = featureDistance(featureReduction)
% featureReduction is a feature array of reduced features

% a simple way of finding differences between syllables based on feature
% alone
% we note that duration is NOT automatically a feature, or it is only one
% of many
    flds = fieldnames(featureReduction);
    featTable = zeros(numel(featureReduction), numel(flds));
    for ii = 1:numel(flds)
        featTable(:,ii) = [featureReduction.(flds{ii})];
    end
    
    % calculating distances with standardized euclidean is a
    % must as we don't standardize our features
    % but we can't use malahanobis b/c we might have NaNs and/or singular
    % values
    corrMatrix = squareform(pdist(featTable,'seuclidean'));
    
    % normalize to [0,1] with linear scaling
    corrMatrix = corrMatrix / max(corrMatrix(:));
end