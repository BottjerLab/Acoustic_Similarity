function [syllables, features, featureRedux] = parseRegionsIntoSyllables(songStruct, regions, params, varargin)
%PARSEREGIONSINTOSYLLABLES divides multiple subregions within a region
%
%   syllables = parseRegionsIntoSyllables(songStruct, regions) finds all
%   the subregions of significantly larger amplitude than their neighbors, 
%   within a previously collected set of regions. This is a batch script
%   which finds a noise filter if necessary and then divides each region
%   one at a time.  
%   
%   [syllables, features] = parseRegionsIntoSyllables(songStruct, regions)
if nargin < 3 || isempty(params),
    params = defaultParams;
end
params = processArgs(params, varargin{:});
fs = 1/songStruct.interval;
params.fs = fs;
%% find a noise region
if params.doFilterNoise && ~isfield(params,'noiseFilter')
    candidateNoise = autodetectNoise(songStruct, regions, params);
    if params.verbose
        fprintf('Done.\n');
        fprintf('Noise segment found, getting noise profile ... ');
    end
    
    % get noise profile
    params.noiseFilter = noiseAnalysis(songStruct,candidateNoise);
    if params.verbose
        fprintf('Done.\n');
    end
end

%% run syllable finding in batch
if params.verbose
    fprintf('Separating syllables from each region...\n');
end
syllables = initEvents;
features = [];
featureRedux = [];
if ~params.plot, progressbar('Breaking down regions...'); end;

nRegions = numel(regions);
for ii = 1:nRegions
    if ~params.plot, progressbar(ii/nRegions); end;
    
    if mod(ii,80) == 0, fprintf('\n'); end; % memory report
    
    clip = getClipAndProcess(songStruct, regions(ii), params);
    rolledRegion = addPrePost(regions(ii),params);
    
    if ~params.quiet, playSound(clip,fs); end

    specType = params.parseSpecType;
    %% spectral analysis
    params.(specType).fs = fs;
    fullSpectrum = getMTSpectrumStats(clip, params.(specType));

    %% segmentation
    protoSylls = segmentSyllables(fullSpectrum, params, 'verbose', false,...
       'optGraphs',{'waveform','deriv','totalPower','wienerEntropy','fracDiffPower'});  
 
%    protoSylls = findPossibleSounds(fullSpectrum, params, 'verbose', false, ...
%        'powerThresh', params.syllable.minPower, ...
%        'riseThresh',params.syllable.minPower * fs / 10);
    
    newSylls = mergeGaps(protoSylls, params.syllable.comboLength/1000);
    nMerged = numel(protoSylls) - numel(newSylls);
    if isempty(newSylls), continue; end;
    
    isTooShort = [newSylls.stop] - [newSylls.start] < params.syllable.minLength/1000;
    newSylls = newSylls(~isTooShort);    
    if params.verbose
    fprintf('Merging %d regions to adjacent ones, %d shorter than %0.2fs...\n', ...
        nMerged, sum(isTooShort), params.syllable.minLength/1000);
    end
    if params.plot, 
        [hax, optGs] = plotAllFigures(fullSpectrum, newSylls, params);
            
        % plot threshold for syllable detection
        axes(hax(strcmp('totalPower', optGs)));
        hold on;
        semilogy(xlim, params.syllable.minPower * [1 1],'b--');
        semilogy(xlim, params.syllable.borderRise * [1 1],'y--');
        
        hold off;
        axes(hax(1));
    end
    
    if params.verbose
        fprintf('Found %d internal syllable%s...\n',numel(newSylls),'s'*(numel(newSylls)~=1));
    end
             
    if params.plot, 
        title(sprintf('Syllable Breakdown (%d/%d) [%s/%s]',ii, nRegions,...
                timeToString(regions(ii).start), ...
                timeToString(regions(ii).stop))); 
        drawnow;
    end;
    
    if ~params.quiet, playSound(clip, fs); end    
    if params.pause, pause; end
    
    % just skip the post processing if there's no syllables found
    if isempty(newSylls), continue; end;
        
    % make sure to adjust time stamps whenever you get regions on
    % subsamples
    absSylls = adjustTimeStamps(newSylls, rolledRegion.start);
    [absSylls(:).type] = deal(2); % index for syllable
    
    syllables = [syllables; absSylls];

    if nargout >= 2
        fieldsToRemove = intersect(fieldnames(fullSpectrum),{'psd','spectrum','freqs','deriv'});
        % while we have the isolated clip, extract features from the clip
        newFeatureRedux = extractFeatures(fullSpectrum, newSylls);
        features = [features rmfield(fullSpectrum,fieldsToRemove)];
    end
    if nargout >= 3
        featureRedux = [featureRedux newFeatureRedux];
    end
    % any filtering should go after this, now that all the features are
    % computed
end