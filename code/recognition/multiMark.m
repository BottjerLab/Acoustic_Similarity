function [marks, regions] = multiMark(songStruct, regions, params, varargin)

if nargin < 3
    params = defaultParams;
end
params = processArgs(params,varargin{:});

if ~isEvent(regions)
    error('markRegions:badArgs', 'regions should be an event structure');
end

if ~isempty(songStruct)
    fs = 1/songStruct.interval;
else
    [~,fs] = getClip(regions(1), songStruct);
end

N = '';
while ~isnumeric(N) || ~(N < 10 && N > 0)
    N = input('Number of different marks? ');
end

labels = cell(N,1);
for ii = 1:N
    if (ii == 1)
        labels{ii} = input('Label #1 (default)? ' ,'s');
    else
        labels{ii} = input(sprintf('Label #%d? ',ii),'s');
    end
end

ii = 1;
marks = NaN(1,numel(regions));
thisParams = params;
stepRoll = 25; % ms

helpAsked = false;
while ii <= numel(regions)
    % get clip
    clip = getClip(addPrePost(regions(ii), thisParams), songStruct);
    
    % get player object
    player = playSound(clip, fs, false);
    
    % request and parse input
    timeInfo = sprintf('%0.2fs', numel(clip)/fs);
    
    helpStr = '? for help';
    fullHelpStr = 'm to menu, u to undo, r to replay, j/h to +-pre roll, k/l to +/-post roll, g to jump, ? for help, ! to abort';
    if helpAsked
        helpStr = fullHelpStr;
    end
    
    prompt = sprintf(['Mark #%d/%d [%s] (1-%d to mark, %s, default is 1)? '], ...
        ii, numel(regions), timeInfo, N, helpStr);
    char = input(prompt,'s');
    if isempty(char),
        char = '1';
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
        case 'm'
            fprintf('Label menu: \n');
            for jj = 1:N
                fprintf('\tLabel #%d: [%s]\n', jj, labels{jj});
            end
            continue;
            % case 'y'
            %     marks(ii) = true;
        case '?'
            helpAsked = true;
            continue;
        case '!'
            fprintf('\nExiting...\n');
            break;
        case 'j'
            thisParams.preroll = thisParams.preroll + stepRoll; continue;
        case 'h'
            thisParams.preroll = max(thisParams.preroll - stepRoll,0); continue;
        case 'k'
            thisParams.postroll = thisParams.postroll + stepRoll; continue;
        case 'l'
            thisParams.postroll = max(thisParams.postroll - stepRoll,0); continue;
        case 'g'
            val = floor(input('Which clip to jump to? '));
            assert(val >= 1 && val <= numel(regions));
            ii = val;
            continue;
        case {'1','2','3','4','5','6','7','8','9'}
            marks(ii) = str2num(lower(char));
            regions(ii).type = labels{marks(ii)};
    end
    ii = ii + 1;
    thisParams = params;
end

end