function [ccmean,ccmax] = syllableAllCross(songStruct, syllables, params, varargin)
% currently returns a full matrix
%% parameter handling
if nargin < 3; params = defaultParams; end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;
params.fs = fs;
nSyll = numel(syllables);

%% get spectrograms
progressbar('Getting spectra for syllables...');
for ii = 1:nSyll
    progressbar(ii/nSyll);
    %fprintf('Getting spectrum for syllable %d...\n',ii);
    [foo,spectra{ii}] = cleanClipFeats(syllables(ii));
end

%% get co-correlation matrix
progressbar('Outer index', 'Inner Index');
ccmean = zeros(nSyll);
ccmax  = zeros(nSyll);
for ii = 1:nSyll
    for jj = ii+1:nSyll
        %fprintf('Examining pair (%d,%d)...\n',ii,jj);
        %ccResults = similarityScan(spectra{ii},spectra{jj}, params);
        ccResults = standardDistance(spectra{ii},spectra{jj}, params);
        ccmean(ii,jj) = mean(ccResults);
        ccmean(jj,ii) = ccmean(ii,jj); % just to fill it
        ccmax(ii,jj)  = max(ccResults);
        ccmax(jj,ii)  = ccmax(ii,jj);
        progressbar([],(jj - ii)/(nSyll-ii));
    end
    progressbar(ii/nSyll,[]);
    ccmean(ii,ii) = 1.0; % just to fill it
    ccmax(ii,ii) = 1.0;
end

%% make reductions

% ordinalize
%cclist = sort([cc{:}]);
%ccOrdMean = arrayfun(@(x) find(x <= cclist,1), ccmean);
%ccOrdMax = arrayfun(@(x) find(x <= cclist,1), ccmax);

% use anchor syllables
%% display
% subplot(2,1,1);
% imagesc(ccmean)
% colorbar
% title('mean')
% subplot(2,1,2);
% imagesc(ccmax)
% colorbar
% title('max')
% subplot(2,2,3);
% imagesc(ccOrdMean)
% colorbar
% title('mean ordinalize')
% subplot(2,2,4);
% imagesc(ccOrdMax)
% colorbar
% title('max ordinalize')

    function [clip,spec] = cleanClipFeats(region)
        clip = getClip(addPrePost(region,params), songStruct);
        
        if ~params.quiet, playSound(clip,fs); end
        
        % high pass the result
        params.fs = fs;
        clip = highPassSample(clip,params);
        
        %% spectral analysis
        params.fine.fs = fs;
        spec = getMTSpectrumStats(clip, params.fine);
        fieldsToRemove = intersect(fieldnames(spec),{'psd','spectrum','freqs','deriv'});    
        spec = rmfield(spec,fieldsToRemove);
    end
end