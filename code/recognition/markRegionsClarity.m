function marks = markRegionsClarity(songStruct, regions, params, varargin)

if nargin < 3
    params = defaultParams;
end
params = processArgs(params,varargin{:});

if ~isEvent(regions)
    error('markRegions:badArgs', 'regions should be an event structure');
end

speedUp = 0;
fsOriginal = 1/songStruct.interval;

ii = 1;
marks = false(1,numel(regions));
thisParams = params;
stepRoll = 25; % ms

helpAsked = false;
speedFac = 1.2;
while ii <= numel(regions)
    % get clip
    clip = getClipAndProcess(songStruct, addPrePost(regions(ii), thisParams), thisParams);
    
    % get player object
    fs = fsOriginal * speedFac ^ speedUp;
    player = playSound(clip, fs, false);
    
    % request and parse input
    timeInfo = sprintf('%0.2fs', numel(clip)/fsOriginal);
    if speedUp ~= 0
        timeInfo = sprintf('%s, %0.1fx', timeInfo, speedFac ^ speedUp);
    end
    
    helpStr = '? for help';
    fullHelpStr = 'y to mark, u to undo, r to replay, j/h to +/-pre roll, k/l to +/-post roll, q/w to speed up/slow down, ? for help';
    if helpAsked
        helpStr = fullHelpStr;
    end
    prompt = sprintf(['Mark #%d/%d [%s] (%s, default is no mark) ? '],ii,numel(regions), timeInfo, helpStr);
    char = input(prompt,'s');
    if isempty(char)
        char = 'n';
    end
    if length(char) > 1
        char = char(1);
    end
    
    helpAsked = false;
    % kill the soundplayer to allow shortcutting behavior
    stop(player);
    
    switch lower(char)
        case 'u'
            if ii > 1, ii = ii - 1; end; thisParams = params; continue;
        case 'r'
            continue;
        case '?'
            helpAsked = true;
            continue;
        case 'y'
            marks(ii) = true;
        case 'j'
            thisParams.preroll = thisParams.preroll + stepRoll; continue;
        case 'h'
            thisParams.preroll = max(thisParams.preroll - stepRoll,0); continue;
        case 'k'
            thisParams.postroll = thisParams.postroll + stepRoll; continue;
        case 'l'
            thisParams.postroll = max(thisParams.postroll - stepRoll,0); continue;
        case 'q'
            speedUp = speedUp + 1; continue;
        case 'w'
            speedUp = speedUp - 1; continue;
    end
    
    ii = ii + 1;
    thisParams = params;
end

end