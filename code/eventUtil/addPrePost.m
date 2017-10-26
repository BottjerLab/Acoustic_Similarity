function regions = addPrePost(regions, params, varargin)

if nargin < 2 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});

if ~isEvent(regions), error('addPrePost:badArgs','Not an event'); end;
% fs = region.idxStart/region.start;
% these are similar to moving a region, but currently they do not accept
% negative prerolls
% assumes events has valid fs

for ii = 1:numel(regions)
    fs = regions(ii).idxStop / regions(ii).stop; 

    regions(ii).start = max(0, regions(ii).start - max(params.preroll / 1000,0));
    regions(ii).idxStart = max(1, regions(ii).idxStart - floor(max(fs * params.preroll / 1000, 0)));

    % note, this doesn't check for the end of the sample, have getClip do this     
    regions(ii).stop = regions(ii).stop + max(params.postroll / 1000,0);
    regions(ii).idxStop = regions(ii).idxStop + ceil(max(fs * params.postroll / 1000,0));
end