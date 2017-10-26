function hline = plotWaveform(clip, fs)
%PLOTWAVEFORM  plots waveform in the window, with downsampling if necessary
len = length(clip);
xx = 1:len;

% Plots downsampled data. This has much faster response time. This is
% based on dsplot by Jiro Doke in file exchange.

% Find outliers
filterWidth = min([50, ceil(length(xx)/10)]); % max window size of 50
a  = clip - ...
    filter(ones(filterWidth,1)/filterWidth, 1, clip);
this.iOutliers = ...
    find(abs(a - repmat(mean(a), size(a, 1), 1)) > ...
    repmat(4 * std(a), size(a, 1), 1));

range = floor([0 len]);

% Despite the spirit of this script (which is to keep memory down for
% large (1M #) arrays,
% We keep a large number of points in order to make the playback simpler -
% since the data needs to be stored anyway for playback in some
% applications (editEventLabels), we might as well keep them all.
% TODO: make this an option in the future
numPoints = 1e6;

if length(xx) > numPoints
    idx = round(linspace(xx(1), xx(end), numPoints))';
    idxi = unique([idx; this.iOutliers]);
    hline = plot(idxi/fs,clip(idxi));
else
    hline = plot(xx/fs,clip(xx));
end

%%%
maxNTicks = 25;
interval = max(0.05,numel(clip)/fs/maxNTicks);
interval = round(interval * 100)/100; % fix to two decimal points
set(gca,'XTick',0:interval:numel(clip)/fs);
xlim([0 len]/fs);
ylim([-1 1]);

end

% Plots downsampled data. This has much faster response time. This is
% based on dsplot by Jiro Doke in file exchange.
% TODO Add data cursor callbacks
function dsplot(clip)
len = length(clip);
x = 1:len;
% Find outliers
filterWidth = min([50, ceil(length(x)/10)]); % max window size of 50
a  = clip - ...
    filter(ones(filterWidth,1)/filterWidth, 1, clip);
[this.iOutliers, this.jOutliers] = ...
    find(abs(a - repmat(mean(a), size(a, 1), 1)) > ...
    repmat(4 * std(a), size(a, 1), 1));

range = floor([0 len]);
numPoints = 50000;
%len = length(this.ChannelHandles);
x = (1:size(clip, 1))';
id = find(x >= range(1) & x <= range(2));
if length(id) > numPoints/len
    idx = round(linspace(id(1), id(end), round(numPoints/len)))';
    ol = this.iOutliers > id(1) & this.iOutliers < id(end);
    for i=1:len
        idxi = unique([idx; this.iOutliers(ol & this.jOutliers == i)]);
        set(this.ChannelHandles(i), 'XData', idxi/fs, ...
            'YData', clip(idxi));
    end
else
    for i=1:len
        set(this.ChannelHandles(i), 'XData', id/fs, ...
            'YData', clip(id));
    end
end
end