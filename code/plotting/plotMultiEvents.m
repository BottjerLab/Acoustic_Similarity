function plotMultiEvents(songStruct, lims, varargin)
    fs = 1/songStruct.interval;
    totalTime = length(songStruct.values) / fs;
    colors = {'r.-','g.-','b.-','k'};
    for ii = 1:numel(varargin)
        evs = varargin{ii};
        xs = [[evs.start]; [evs.stop]; nan(1,numel(evs))];
        ys = [ones(2,numel(evs)) * ii ; nan(1,numel(evs))];
        plot(xs(:),ys(:),colors{mod(ii-1,numel(colors))+1},'MarkerSize', 8,'LineWidth', 2   );
        hold on;
    end
    xlim([0 totalTime]);
    ylim([0.5 0.5 + numel(varargin)])
    hold off;
end