function [vocString, clusterIdxs] = createAlphabet(syllables, distMatrix, songStruct, params, varargin)
%CREATEALPHABET given distance matrix, cluster into groups
% 
if nargin < 4 || isempty(params)
    params = defaultParams;
end
params = processArgs(params,varargin{:});

fs = 1/songStruct.interval;
% cluster the syllables by distance - assumes triangular matrix
%boutCorrVector = squareform(corrMatrix + corrMatrix' ...
%    - 2 * diag(diag(corrMatrix)),'tovector');
% dendTree = clusterAll(vocalizations, boutCorrMatrix, songStruct);
%dendTree = linkage(1 - boutCorrVector, 'complete'); % TODO: try top-down instead of agglomerative?

% for correlation matrices
% dendTree = linkage(1 - corrMatrix, params.clusterMethod); 
% for distance matrices
dendTree = linkage(distMatrix, params.clusterMethod);
%% create clusters

clusterIdxs = cluster(dendTree, 'maxClust', params.nClusters)';

%% reassign numbers to cluster letters - don't do this for now
alphaFreq = histc(clusterIdxs,1:params.nClusters);
% translate to string
vocString = numToAlpha(clusterIdxs);

freqs = hist(clusterIdxs, 1:params.nClusters)
if params.plot
    asciied = double(vocString);
    hist(asciied,min(asciied):1:max(asciied))
    xlim([min(asciied)-0.5 max(asciied)+0.5])
    set(gca,'XTick',min(asciied):1:max(asciied))
    set(gca,'XTickLabels',char(min(asciied):1:max(asciied))');
    title('Distribution of labels');
end
% sample the alphabet
if params.playsample
    
    for ii = 1:params.nClusters
        % pick at most 5 example syllables and more if we have more
        nSamples = max(floor(sqrt(alphaFreq(ii))),min(5,alphaFreq(ii)));
        fprintf('Playing syllable %s (# = %d)...\n',numToAlpha(ii), freqs(ii));
        sampleIdxs = find(clusterIdxs==ii,nSamples);
        for jj = 1:nSamples
            playSound(getClip(syllables(sampleIdxs(jj)), songStruct),fs,true);
            pause(0.25);
        end
        beep
        pause(0.25);
    end
end

end

function str = numToAlpha(vec)
    vec(vec>= 27 & vec <=52) = vec(vec>= 27 & vec <=52) + double('A') - double('a') - 26;
    str = char(double('a' - 1 + vec));
end

function vec = alphaToNum(str)
    vec = str - 'A' + 1;
end