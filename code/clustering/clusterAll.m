function links = clusterAll(syllables, corrMatrix, songStruct, params, varargin)
% function links = clusterAll(syllables, corrMatrix, songStruct, params, varargin)
%
% simple greedy agglomerative (bottom-up clustering)
% takes correlation matrix, and clusters syllables
% by similarity as defined in the corrMatrix
% corrMatrix is a block-diagonal matrix, each block represents a
% bout/sub-part of song

if nargin < 5 || isempty(params)
    params = defaultParams;
end
params = processArgs(params,varargin{:});

nSylls = numel(syllables);

% closest objects
subLinks{ii} = zeros(nSylls - 1,3);
parent = 1:nSylls;
nInTree = ones(1,nSylls);
distances = 1 - subCorr;

% this clustering is very similar to the linkage function, but
% simple to write
for jj = 1:nSylls - 1
    % clI < clJ b/c it's upper triangular
    [minval, idx] = nanmin(distances(:));
    [clI, clJ] = ind2sub(size(distances), idx);
    
    subLinks{ii}(jj,1) = parent(clI);
    subLinks{ii}(jj,2) = parent(clJ);
    subLinks{ii}(jj,3) = minval;
    
    linkageType = 'unweighted';
    % update weights according to link function...
    % NaN the lower indices
    for kk = 1:nSubIdxs
        % for every other idx,
        % take the maximum of the two previous groups
        if clJ ~= kk
            distI = distances(min(clI,kk),max(clI,kk));
            distJ = distances(min(clJ,kk),max(clJ,kk));
            switch linkageType
                case 'single'
                    distances(min(clJ,kk),max(clJ,kk)) = min(distI,distJ);
                case 'complete'
                    distances(min(clJ,kk),max(clJ,kk)) = max(distI,distJ);
                case 'weighted'
                    distances(min(clJ,kk),max(clJ,kk)) = 1/2 * (distI + distJ);
                case 'unweighted'
                    distances(min(clJ,kk),max(clJ,kk)) = ...
                        (distI * nInTree(clI) + distJ * nInTree(clJ)) / (nInTree(clI) + nInTree(clJ));
            end
        end
        distances(min(clI,kk),max(clI,kk)) = NaN;
    end
    
    % update parentage
    parent(clI) = nSubIdxs + jj;
    parent(clJ) = nSubIdxs + jj;
    nInTree(clI) = nInTree(clI) + 1;
    nInTree(clJ) = nInTree(clJ) + 1;
end


end
