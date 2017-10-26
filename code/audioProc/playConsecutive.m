function playConsecutive(songStruct, events, params, varargin)
    if nargin < 3
        params = defaultParams;
    end
    if nargin > 3
        params = processArgs(params, varargin{:});
    end
    
    fs = 1/songStruct.interval;
    for ii = 1:numel(events)
        playSound(getClipAndProcess(songStruct, events(ii), params),fs);
    end
end