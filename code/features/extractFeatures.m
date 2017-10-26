function features = extractFeatures(spectrum, regions)
%EXTRACTFEATURES reduce a feature (TODO)
%
%    features = extractFeatures(spectrum, regions)
% summarize the features out of a spectrum from a clip
% regions must be relative to the clip start from which the spectrum was
% taken
% if spectrum doesn't have any data, returns a feature vector of all 0's, except for length.

%   Parameter options:
%      highPassFreq: The frequency above which the sound is high-passed:
%      fine: the spectrum subparameters over which to collect the spectrogram
%      plot: [t/f] controls the plotting output
%      pause: [t/f] controls pausing during plotting

%% set defaults
if nargin < 2 && numel(spectrum.times) > 0
    regions.start = spectrum.times(1);
    regions.stop = spectrum.times(end);
elseif nargin < 2, % just the one clip
    regions.start = 0; regions.stop = 0;
end

%% find all features that can be summarized
fields = whichFeatures(spectrum);

%% create empty structure with same fields
% these must all be reductions, vector to scalar (i.e. R^n->R)
suffixes = {'nanmean', 'nanstd', 'max', 'min', ...
    'attack','sustain','release',... % custom functions
    'lintrend', 'peak', 'trough'}; 

% define the functions for attack, sustain and release (note: we require
% no NaNs here)
    function ret =   attack(vec), ret = vec(ceil(end * 0.1)); end %we must change this average of sample 1-3%
    function ret =  sustain(vec), ret = vec(ceil(end * 0.5)); end %we must change this average of sample 49-51%   
    function ret =  release(vec), ret = vec(ceil(end * 0.9)); end %we must change this average of sample 97-100% 
    function ret = lintrend(vec), ret = corr((1:length(vec))', vec'); if isnan(ret), ret = 0; end; end
    function ret =     peak(vec), [~, ret] = max(vec); end
    function ret =   trough(vec), [~, ret] = min(vec); end

% a trick to tell MATLAB these are nested functions
fxnNames = suffixes; % all reductions
fxnNames(5:end) = strcat('extractFeatures/', fxnNames(5:end)); 

% warning; may not be stable past before R2012b
% this converts the names into functions, for brevity
fhs = cellfun(@str2func, fxnNames, 'UniformOutput',false);

%TODO: check NaN behavior with the boosting classification framework

statFields = ['length'];
for ii = 1:numel(fields)
    statFields = [statFields strcat(fields(ii), strcat('_', suffixes))];
end
features = initEmptyStructArray(statFields);

for ii = 1:numel(regions)
    if numel(spectrum.times) > 0
        st = find(regions(ii).start <= spectrum.times,1);
        en = min(find(regions(ii).stop  <= spectrum.times,1));
        if isempty(en), en = numel(spectrum.times); end;

        for jj = 1:numel(fields) 
            % grab each field's data
            featData = spectrum.(fields{jj})(st:en);
            if numel(featData) == 0
                keyboard
                error('extractFeatures:shortSpectrum','Empty feature [%s] data due to short spectrum', fields{jj});
            end
            
            for kk = 1:numel(suffixes)  % apply each reduction to the data
                features(ii).([fields{jj} '_' suffixes{kk}]) = fhs{kk}(featData);
            end
        end
    else
        for jj = 1:numel(fields)
            for kk = 1:numel(suffixes)
                features(ii).([fields{jj} '_' suffixes{kk}]) = 0;
            end
        end
    end
    features(ii).length = regions(ii).stop - regions(ii).start;
end
end