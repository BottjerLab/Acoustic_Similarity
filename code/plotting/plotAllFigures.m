function [hax, optGraphs, hfig] = plotAllFigures(spec, regions, params, varargin)
% PLOTALLFIGURES main plotting function
%   plotAllFigures(spec) plots the waveform in parallel with audio features
  
% handling bad arguments 
% (1) empty spectrogram
if isempty(spec)
    hax = 0; 
    optGraphs = params.optGraphs;
    hfig = 0; 
    return; 
end % returns axis handles

if nargin < 2, regions = []; end; 
if nargin < 3 || isempty(params), params = defaultParams; end;

% plot only the graphs that you want and that have information provided by
% the spectrum data structure
params = processArgs(params, varargin{:});
optGraphs = params.optGraphs;
optGraphs = optGraphs(isfield(spec, optGraphs) | strcmp('spectrogram',optGraphs) | ...
    strcmp('fracDiffPower',optGraphs));

% find the number of plots
nPlots = numel(optGraphs);
hax = zeros(1,nPlots);

% labels for events structure
toggleLabels = true;

cP = 0;
for ii = 1:nPlots
    % plot in next window
    cP = cP + 1;
    hax(ii) = subplot(nPlots,1,cP);
    switch(optGraphs{ii})
        case 'spectrogram'
            plotSpectrogram(spec);
            %     if isfield(spec,'fundamentalFreq')
            %         hold on; plotPitch('fundamentalFreq','r-'); hold off;
            %     end;
            %     if isfield(spec,'centerFreq')
            %         hold on; plotPitch('centerFreq','y-'); hold off;
            %     end;
            for ii = 1:numel(regions), 
                regionsMark(ii).start = regions(ii).start * numel(spec.times) / max(spec.times);
                regionsMark(ii).stop  = regions(ii).stop  * numel(spec.times) / max(spec.times);
            end
            %plotMarks(regions);
        case 'deriv'
%            if any(strcmp(optGraphs,'spectrogram'))
%                freezeColors(hh); %seems not to work in 2011b =/
%            end
            
            % todo: it would be a very nice detail to change contrast based
            % on the scroll wheel, but that would mean carrying large data
            % amounts in the userData memory
            plotDerivGram(spec,params);
            %{
            if isfield(spec,'fundamentalFreq')
                hold on; plotPitch('fundamentalFreq','r-.', 'LineWidth',0.5); hold off;
            end;
            if isfield(spec,'centerFreq')
                hold on; plotPitch('centerFreq','y-.', 'LineWidth',0.5); hold off;
            end;
%}
            switch params.showDgramRegionStyle
                case 'lines' % non-fancy non-interactive method
                    plotMarks(regions);
                    plotLabels;
                case 'fancy' % fancy interactive method, slows down interface some and
                    set(gcf,'Renderer','OpenGL');
                    plotAreaMarks(regions, [0 0.8 0.8 0.15]); % requires openGL
                    plotLabels;
                    % note: have to use opengl software to render, which is
                    % slower
            end
        case 'totalPower'
            plotRMS;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'fracDiffPower'
            plotFracDiffPower;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'wienerEntropy'
            plotEntropy;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'FM'
            plotFM;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'AM'
            plotAM;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'rawAM'
            plotRawAM;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'pitchGoodness'
            plotPitchGoodness;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'aperiodicity'
            hold off; plotFracPeriodicPower;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'waveform'
            if ~isfield(params,'fs')
                error('plotAllFigures:undefined','sampling rate not defined, pass parameters');
            end;
            hline = plotWaveform(spec.waveform,params.fs);
            
            % magic user interface happens here
            set(hline,'UserData',params.fs); % give playback clue
            xlim([min(spec.times) max(spec.times)]);
            hp = plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            
            
            % magic playback UI control happens here
            set(gcf,'WindowButtonDownFcn',{@buttonDownFcn, [hax(ii) hp]});
            toggleLabels = false;
        case 'centerFreq'
            hold off; plotPitch('centerFreq');            
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'fundamentalFreq'
            hold off; plotPitch('fundamentalFreq');
            title('Fundamental frequency');
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'mTD'
            hold off; plotPartialDerivs;
            plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        case 'mfcc'
            hold off; plotMFCC;
            %plotAreaMarks(regions,[],params.showLabels && toggleLabels);
            toggleLabels = false;
        otherwise
            error('plotAllFigures:UnspecifiedPlotting','Don''t know how to plot %s', optGraphs{ii});
    end
end

subplot(nPlots,1,1); % for correct title positioning

% allow for linked scanning/zooming on segments
set(gcf,'UserData', struct('hlink', linkprop(hax, 'XLim')));
%drawnow 
hfig = gcf;

% individual plotting for different acoustic features
    function plotRMS
        % technical note:
        % openGL rendering does not support log scaling, which causes 
        % problems with highlighting 
        
        % plot unsmoothed power 
        semilogy(spec.times, spec.totalPower, 'm-');

        % plot first derivative 
        smoothedPower = smoothSignal(spec.totalPower,17);
        dt = spec.times(2) - spec.times(1);
        dPdT = [diff(smoothedPower) / dt eps];
        posTimes = spec.times; posTimes(dPdT<=0) = NaN;
        negTimes = spec.times; negTimes(dPdT>=0) = NaN;
        
        % plot negative and positive as diff signs
        hold on;
        semilogy(posTimes, dPdT,'g-')
        semilogy(negTimes, -dPdT,'r-')
        hold off;
        
        % set plotting window
        set(gca, 'YLimMode', 'manual')
        minPower = min(spec.totalPower);
        ylim([minPower/2,1e-1])
        set(gca,'YTick',10.^(ceil(log10(minPower)):-1));
        xlim([min(spec.times) max(spec.times)]);
        
        ylabel('log RMS Power')
        title('Total Power (Volume)')
    end

    function plotFracDiffPower
        
        % plot first derivative 
        smoothedPower = smoothSignal(spec.totalPower,17);
        dt = spec.times(2) - spec.times(1);
        dPdT = [diff(smoothedPower) / dt eps];
        dfPdT = dPdT./spec.totalPower;
        posTimes = spec.times; posTimes(dPdT<=0) = NaN;
        negTimes = spec.times; negTimes(dPdT>=0) = NaN;

        semilogy(posTimes, dfPdT,'g-', negTimes, -dfPdT,'r-')
        hold off;
        
        % set plotting window
        set(gca, 'YLimMode', 'manual')
        xlim([min(spec.times) max(spec.times)]);
        %minPower = min(spec.totalPower);
        ylim([min(abs(dfPdT))/3,max(abs(dfPdT)) * 3])
        %set(gca,'YTick',10.^(ceil(log10(minPower)):-1));

    end
    function plotEntropy
        plot(spec.times, spec.wienerEntropy, 'r-');
        xlim([min(spec.times) max(spec.times)]);
        ylabel('Entropy')
        title('Wiener entropy (high = noise, low = tone, mid = stack)')
    end
    function plotFM
        %semilogy(spec.times, spec.mTD, 'r-',spec.times, spec.mFD,'b-');
        
        plot(spec.times, spec.FM, 'r-');
        xlim([min(spec.times) max(spec.times)]);
        %ylim([0 90]);
        ylabel('FM (deg)')
        title('Frequency modulation')
    end
    function plotAM
        plot(spec.times, spec.AM, 'r-');
        xlim([min(spec.times) max(spec.times)]);
        %ylim([0 90]);
        ylabel('AM (1/ms)')
        title('Amplitude modulation')
    end
    function plotRawAM
        % plot first derivative 
        posTimes = spec.times; posTimes(spec.rawAM<=0) = NaN;
        negTimes = spec.times; negTimes(spec.rawAM>=0) = NaN;
        
        % plot negative and positive as diff signs
        semilogy(posTimes, spec.rawAM,'g-')
        hold on;
        semilogy(negTimes, -spec.rawAM,'r-')
        hold off;
        
        xlim([min(spec.times) max(spec.times)]);
        ylabel('raw AM (dB/ms)')
        ylim([1e-5 1]);
        title('Unnormalized amplitude modulation')
    end
    function plotPitchGoodness
        semilogy(spec.times, spec.pitchGoodness, 'r-');
        % plot(spec.times, log10(spec.pitchGoodness), 'r-');
 
        xlim([min(spec.times) max(spec.times)]);
        ylabel('log Cepstral strength');
        title('Goodness of pitch');
    end
    function plotPitch(field, colStyle, varargin)
        if nargin < 2, colStyle = 'r-'; end;
        if nargin < 1, field = 'fundamentalFreq'; end;

        plot(spec.times, spec.(field), colStyle,varargin{:});
        xlim([min(spec.times) max(spec.times)]);
        ylim([min(spec.freqs) max(spec.freqs)]);
        ylabel('freq (Hz)');
%        title(field);
    end
    function plotFracPeriodicPower
        plot(spec.times, spec.aperiodicity, 'r-');
        xlim([min(spec.times) max(spec.times)]);
        ylim([0 1]);
        ylabel('aper pow frac');
        title('Fraction of Periodic power')
    end
    function plotPartialDerivs
        semilogy(spec.times, spec.mTD, 'r-', spec.times, spec.mFD, 'g-');
        %plot(spec.times, log10(spec.mTD), 'r-', ...
        %     spec.times, log10(spec.mFD), 'g-');
        xlim([min(spec.times) max(spec.times)]);
        ylabel('partial deriv amp (dB/msec)');
    end
    function plotMFCC
        nColors=256;
        % the first row of the MFCC is the volume
        nBands = size(spec.mfcc, 1) - 1;
        hs = imagesc([min(spec.times) max(spec.times)],...
             [2 nBands],... 
             spec.mfcc(2:end,:));
        xlabel('Time (s)');
        ylabel('mel-freq cepstral coeff');
        set(gca,'YDir','normal');
        colormap(gray(nColors));
    end
    function plotLabels
        for kk = 1:numel(regions)
            if ~isempty(regions(kk).type)
                xpos = (regions(kk).start + regions(kk).stop)/2;
                createTextLabel(xpos, regions(kk).type);
            end
        end
    end

end    


function buttonDownFcn(gcbo, eventdata, handles)
hax = handles(1); hp = handles(2);
currPt = get(hax,'CurrentPoint'); currPt = currPt(1,1:2);

[patchClicked, lineClicked, hitWindow] = clickStatus(currPt, handles);
if ~hitWindow, return; end;

% if we are NOT double clicked, get out
if ~strcmp(get(gcbo,'SelectionType'),'open'), return; end;

% we are double clicked: play a sound
% get the clip data

hline = findobj(hax,'Type','line');
hline = hline(end); % it should be the one furthest back
xWave = get(hline,'XData'); yWave = get(hline,'YData');
fs = get(hline,'UserData');
if isempty(patchClicked) % what should we do if a patch is not clicked?
    % play the whole clip
    playSound(yWave, fs, true);
else
    % get the borders
    currXData = get(hp, 'XData');
    patchBorders = currXData(2:3,patchClicked);
    
    % cut the clip at the right points
    clipStart = find(xWave >= patchBorders(1),1);
    clipEnd = find(xWave >= patchBorders(2),1);
    clip = yWave(clipStart:clipEnd);
    
    % play the sound clip, while blocking
    playSound(clip, fs, true);
end
end

function [patchSeld, lineSeld, hitWindow] = clickStatus(currPt, handles)
% returns empties on default
patchSeld = []; lineSeld = [];
hax = handles(1); hp = handles(2);

win([1 3]) = get(hax,'XLim'); win(3) = win(3) - win(1);
win([2 4]) = get(hax,'YLim'); win(4) = win(4) - win(2);
hitWindow = inRect(win, currPt); if ~hitWindow, return, end;

yy = get(hax,'YLim');

if currPt(2) > yy(2) || currPt(2) < yy(1), return; end;

xBounds = get(hp,'XData');
if isempty(xBounds), return; end; % nothing to click
xBounds = xBounds(2:3,:);

patchSeld = find(xBounds(1,:) <= currPt(1) & xBounds(2,:) >= currPt(1));

if numel(patchSeld) > 1
    % find the one which is closer
    distsToCursor = min(xBounds(1,patchSeld) - currPt(1));
    [~,closest] = min(distsToCursor);
    patchSeld = patchSeld(closest);
end
end

function foo = inRect(win, pt)
foo = win(1) <= pt(1) && pt(1) < win(1) + win(3) && ...
    win(2) <= pt(2) && pt(2) < win(2) + win(4);
end

