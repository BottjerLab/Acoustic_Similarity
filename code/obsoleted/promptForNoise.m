function noiseRegion = promptForNoise(songStruct, params, varargin)
% PROMPTFORNOISE interactively marks noise for filtering
%
if ~exist('params') || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});
params.fs = 1 / songStruct.interval;

% this code ensures that we don't cut any songs into two pieces
% we keep an overlap of ~300ms in each clip
wholeSong = getWholeSongAsRegion(songStruct);
windows = splitIntoOverlap(wholeSong, params.Nsplits, params.overlapSplit / 1000);
% LEN_OVERLAP = params.fs * params.overlapSplit / 1000;
% LEN_SEGMENT = (songStruct.length - LEN_OVERLAP) / params.Nsplits + ...
%     LEN_OVERLAP;

noiseRegion = [];
for ii = 1:params.Nsplits
    fprintf('Select noise region...\n');
    noiseRegion = plotAndAdjust(songStruct, [], windows(ii),params);
    if ~isempty(noiseRegion), break; end;
end

end