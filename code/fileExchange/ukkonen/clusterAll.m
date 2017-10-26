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


% make grand dendrogram
% format is 3 x N - 1, first two columns are linked nodes,
% third column is distance
% new nodes are N + i;
subCorr = corrMatrix;
subCorr(tril(true(nSylls),0)) = NaN; % replace the lower half and diagonal with NaNs

% closest objects
links = zeros(nSylls - 1,3);
parent = 1:nSylls;
nInTree = ones(1,nSylls);
distances = 1 - subCorr;

% this clustering is very similar to the linkage function, but
% simple to write
for jj = 1:nSylls - 1
    fprintf('.');
    if mod(jj,80) == 0, fprintf('\n'); end;
    % clI < clJ b/c it's upper triangular
    [minval, idx] = nanmin(distances(:));
    [clI, clJ] = ind2sub(size(distances), idx);
    
    links(jj,1) = parent(clI);
    links(jj,2) = parent(clJ);
    links(jj,3) = minval;
    
    linkageType = 'unweighted';
    % update weights according to link function...
    % NaN the lower indices
    for kk = 1:nSylls
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
    parent(clI) = nSylls + jj;
    parent(clJ) = nSylls + jj;
    nInTree(clI) = nInTree(clI) + 1;
    nInTree(clJ) = nInTree(clJ) + 1;
end

%     %% generate dendrogram plot
%     % 2nd argument in forces all the leaves to be shown
%     if ii == 1 || ~params.plot, hh = figure; end;
%     [foo,hhand,perm] = dendrogram(links, nSylls+1,...
%         'Labels',cellstr(num2str((1:nSylls)')), ...
%         'Orientation','left');
%     title(sprintf('Dendrogram for Bout #%d',ii));
%     if ~params.plot, close(hh); end
%     if nargin >= 4 && params.playsample
%         fs = 1/songStruct.interval;
%         % draw the moving green arrow
%         if params.plot
%             hold on;
%             harrow = plot(xlim * [0.98 0.02]', 0, 'g>','MarkerFaceColor','g');
%         end
%         % play each sample in clustered order
%         for jj = 1:nSylls
%             if params.plot, set(harrow,'YData',jj); end;
%             eve = addPrePost(syllables(idxs{ii}(perm(jj))),...
%                 processArgs(defaultParams,'preRoll',0.02,'postRoll',0.02));
%             clip = getClip(songStruct, eve);
%             playSound(clip, fs);
%             pause(0.8);
%         end;
%     end
%     if params.plot, hold off; end;
%     if params.plot
%         pause(1.0);
%     end

end
