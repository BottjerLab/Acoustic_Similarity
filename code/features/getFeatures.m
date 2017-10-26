function [features,paddedSpecs] = getFeatures(songStruct, regions, params, varargin)
% returns feature and spectra 
% note: can take empty songStruct if syllables are DR (disk read)
if nargin < 3 || isempty(params)
    params = defaultParams;
end
if numel(varargin) > 0
    params = processArgs(params, varargin{:});
end
if ~isempty(songStruct)
    params.fs = 1/songStruct.interval;
else
    % params.fs must be preset if using disk read syllables
    assert(isfield(params,'fs') && ~isnan(params.fs));    
end
%% preprocessing
% filter out noise first

%params = processArgs(params,'preroll', 3, 'postroll', 3);

nRegions = numel(regions);
%progressbar(sprintf('Gathering features (#=%d)',nRegions));
for ii = 1:nRegions
    roi = regions(ii);
	paddedClip = getClipAndProcess(songStruct, roi, params);
       
    %% run spectra analysis 
    params.fine.fs = params.fs;

    
    tempSpec = getMTSpectrumStats(paddedClip, params.fine);
    if isempty(tempSpec), continue; end; 
    % the tempSpec may return empty
    
    % note: take out the harmonicPitch field because it's dependent on the
    % sampling rate

    %% play sample if required
    if params.playsample
        playSound(paddedClip, params.fs, true);
    end
    
    %% get features
    featSpecs = initEmptyStructArray(params.reduceFeatures, 1);
    for jj = 1:numel(params.reduceFeatures)
        featSpecs.(params.reduceFeatures{jj}) = tempSpec.(params.reduceFeatures{jj});
    end
    featSpecs.times = tempSpec.times;
    regionInContext = adjustTimeStamps(roi, -roi.start);
    features(ii) = extractFeatures(featSpecs, regionInContext);
    
    paddedSpecs{ii} = featSpecs;
    if params.plot
        plotAllFigures(paddedSpecs(ii), regionInContext, params);
    end
 %   progressbar(ii/nRegions);
end