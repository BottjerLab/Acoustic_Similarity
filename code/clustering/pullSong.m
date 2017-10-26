function ccmean = pullSong(songStruct, events, syllables, specs, params, varargin)
% NB: candidate for removal, similar to syllableAllCross.m

%% parameter handling
if nargin < 5; params = defaultParams; end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;

% get the max syllable number
maxSylIdx = find(max([events.idxStop]) <= [syllables.idxStop], 1);
ccmean = zeros(maxSylIdx);
%% get spectrograms
for hh = 1:numel(events)
    fprintf('Working on event %d...\n',hh);
    [subSylls, subidxs] = getSubEvents(events(hh), syllables);
    nSub = numel(subSylls);
    if nargin < 4
        spec = [];
        for ii = 1:nSub
            fprintf('Getting spectrum for syllable %d...\n',ii);
            [foo,spec(ii)] = cleanClip(subSylls(ii));
        end
    else
        spec = specs(subidxs);
    end
    cc = cell(nSub);
    %% get cocorrelation matrix
    for ii = 1:nSub
        for jj = ii+1:nSub
            fprintf('Examining pair (%d,%d)...\n',ii,jj);
            cc{ii,jj} = similarityScan(spec(ii),spec(jj));
            cc{jj,ii} = [0.0]; % just to fill it
        end
        cc{ii,ii} = [0.0]; % just to fill it
    end
    
    %% make reductions
    ccmean(subidxs, subidxs) = cellfun(@(x) mean(x),cc);
end
imagesc(ccmean);
% ordinalize
%{
cclist = sort([cc{:}]);
ccOrdMean = arrayfun(@(x) find(x <= cclist,1), ccmean);
ccOrdMax = arrayfun(@(x) find(x <= cclist,1), ccmax);

% use anchor syllables
% in our segmentation, these are long syllables with low mean correlation to
% other syllables but high self correlation.
% these are (a) low entropy - highly peaked and (b) only common with other
% low entropy
entropy = zeros(1,nSub);
ccOrdP = ccOrdMean/numel(cclist);
for ii = 1:nSub
    for jj = 1:ii-1
        ccOrdP(ii,jj) = ccOrdP(jj,ii);
    end
end

for ii = 1:nSub
    for jj = 1:nSub
        if ii == jj, continue; end
        pN = ccOrdP(ii,jj) / sum(ccOrdP(ii,:));
        entropy(ii) = entropy(ii) - pN * log(pN);
    end
end

% see if the least entropic N syllables have high correlation
[foo,sortSyll] = sort(entropy, 1,'descend');
avgAgreement = ones(1,nSub);
for ii = 2:nSub
    avgAgreement(ii) = 0;
    for jj = 1:ii
        for kk = 1:jj-1
            avgAgreement(ii) = avgAgreement(ii) + ccOrdP(sortSyll(jj), sortSyll(kk));
        end
    end
    avgAgreement(ii) = avgAgreement(ii) / (ii * (ii + 1) / 2);
end
%}
%% display


    function [clip,spec] = cleanClip(region)
        clip = getClip(addPrePost(region,params), songStruct);
        
        if ~params.quiet, playSound(clip,fs); end
        
        % high pass the result
        params.fs=fs;
        clip = highPassSample(clip,params);
        
        %% spectral analysis
        params.fine.fs = fs;
        spec = getMTSpectrumStats(clip, params.fine);
    end
end