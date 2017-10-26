function plotEvent(songStruct, ev, regions)
    cl = getClip(ev, songStruct);
    plot( (1:numel(cl)) * songStruct.interval, cl);
    plotMarks(0.95)

function plotMarks(yc)
        hold on;
        if ~isempty(regions)
            plot([regions.start],yc, 'cv',...
                'MarkerFaceColor','c', 'MarkerSize', 5);
            plot([regions.stop],yc, 'mv',...
                'MarkerFaceColor','m', 'MarkerSize', 5);
            
            for jj = 1:numel(regions)
                ymins = ylim;
                if strcmp(get(gca,'YScale'),'log')
                    semilogy(repmat([regions.start],2,1), [ymins(1) yc],'c-');
                    semilogy(repmat([regions.stop],2,1), [ymins(1) yc],'m-');
                else
                    plot(repmat([regions.start],2,1), [ymins(1) yc],'c-');
                    plot(repmat([regions.stop],2,1), [ymins(1) yc],'m-');
                end
            end
            
        end
        hold off;
end
end
