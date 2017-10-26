function [hf, hims] = mosaicDRSpec(DRevents, params, varargin)

hf = [];
if isempty(DRevents)
    warning('mosaicDRspec:noClips', 'No clips input or maxLength is smaller than clip');
return;
end

if nargin < 2 || isempty(params)
    params = defaultParams;
end
params = processArgs(params,varargin{:});

% make sure we only have the clips we want
clLens = [DRevents.stop]-[DRevents.start];
cumLens = cumsum(clLens);
if ~isinf(params.maxMosaicLength) && cumLens(end) > params.maxMosaicLength
    cutoff = find(cumLens < params.maxMosaicLength, 1, 'last');
    cumLens = cumLens(1:cutoff);
    clLens = clLens(1:cutoff);
    DRevents = DRevents(1:cutoff);
end

% extend by any roll
DRevents = addPrePost(DRevents, params);

% make one long clip...
nEv = numel(DRevents);
cl = cell(1, nEv);
fs = zeros(1,nEv);

%params = processArgs(params, 'noroll');  % strict boundaries
for ii = 1:nEv
    [cl{ii}, fs(ii)] = getClipAndProcess([],DRevents(ii), params, 'noroll');
    
    % TEMP: readjust lengths (bugfix until we can fix noiseGate chopping
    % problem)
    clLens(ii) = numel(cl{ii}) / fs(ii);
end
cumLens = cumsum(clLens);

if isempty(fs)
    warning('mosaicDRspec:noClips', 'No clips input or maxLength is smaller than clip');
    return
end
if ~all(fs==fs(1))
    error('mosaicDRSpec:variableSampling', ...
    'Different sampling frequencies in each clip...'); 
end
fs = fs(1);

%%
maxClipLen = params.mosaicRowLength;  %seconds

nLowPlots = floor(cumLens(end)/maxClipLen);
concatCl = cell(1,nLowPlots);
idxRange = 1:find(cumLens < maxClipLen, 1,'last');
params.fine.fs=fs;
refRanges = cell(1,nLowPlots);

ii = 1;
while ~isempty(idxRange) && (idxRange(end) <= nEv || isempty(concatCl))
    % concatenate clip
    concatCl{ii} = vertcat(cl{idxRange});
    refRanges{ii} = idxRange;
    idxRange = (idxRange(end)+1):find(cumLens < maxClipLen + cumLens(idxRange(end)), 1, 'last');
    ii = ii+1;
end

% for display purposes - 1/3 rows is the most proportional for spectrogram
% display
nPlots = max(3,numel(concatCl));


%maxClipLen = max(cellfun(@numel, concatCl));

for ii = 1:min(nPlots, numel(concatCl));
    % calculate and plot spectrogram - todo: more space, less ticks...
    %subplot(nPlots,1,ii);
    % handle all plots to scale - 3rd entry has normalized unit of length
    hh(ii) = subplot('Position',[0 (nPlots-ii)/nPlots (length(concatCl{ii})/fs)/maxClipLen 1/nPlots]);
    spectrum = getMTSpectrumStats(concatCl{ii}, params.fine);
    hims(ii) = plotDerivGram(spectrum,params);
    
    % plot vertical separators
    yy = ylim;
    xpts = cumsum(clLens(refRanges{ii}));
    xpts(end) = []; % don't need the last one b/c it's already at the boundary
   
    hold on;
    for jj = 1:numel(xpts)
        plot(xpts(jj) * [1 1], yy,'w-','LineWidth',2);
    end
    
    % plot a small time scale bar
    if ii == 1
        xx = xlim;
        xSB = xx * [0.95 0.05]';
        fracBar = 0.04;
        xSBEnd = xx * [(0.95 - fracBar) (0.05 + fracBar)]'; 
        ySB = yy * [0.15 0.85]';
        textH = yy * [0.24 0.76]';
        
        SBlen = roundn(maxClipLen * 1000 * fracBar, 10); % 1/25th of the bar in milliseconds, 10 ms length intervals;        
        plot([xSB, xSBEnd], [ySB, ySB], 'g-', 'LineWidth', 2);
        text(xSB, textH, sprintf('%d ms', SBlen),'FontWeight','bold',...
        'HorizontalAlignment','left','VerticalAlignment','bottom');
        hold off;
    end
    % remove axes
    set(gca,'YTickLabel',[],'YTick',[]);
    set(gca,'XTickLabel',[],'XTick',[]);
    
end
% reset for title placement
axes(hh(1))
hf = gcf;
end

function x = roundn(x,y)
% round to the nearest y
x = round(x/y)*y;
end
