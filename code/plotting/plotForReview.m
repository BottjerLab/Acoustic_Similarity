function plotForReview(songStruct, regions, superRegions, params, varargin)
    if nargin < 3, superRegions = 100; end;
    if nargin < 4 || isempty(params), params = defaultParams; end;
    if nargin < 5; varargin = {}; end;
        
    if ~isEvent(superRegions), 
        superRegions = splitIntoOverlap ...
                (getWholeSongAsRegion(songStruct), superRegions);
    end
    plotAndAdjust(songStruct, regions, superRegions, params, varargin{:}, ...
        'readOnly', true, 'showLabels', true);
    % no edited regions
end