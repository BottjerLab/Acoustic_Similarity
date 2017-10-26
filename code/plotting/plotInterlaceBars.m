function plotInterlaceBars(setGray, setRed, bins)
% Gray plotted in gray, Red plotted in red

binTol = 1e-5;
hGray = histc(setGray , bins);
hRed  = histc(setRed, bins);

holdState = ishold;
if all(diff(bins) - (bins(2) - bins(1)) < binTol * (bins(2) - bins(1)))
    bw = bins(2) - bins(1);
    bar(bins, hGray, 0.5, 'FaceColor', [0.5 0.5 0.5]);
    hold on;
    bar(bins+bw/2, hRed, 0.5, 'r');
else
    bw = diff(bins); bw = [bins(1)/2 bw bw(end)];
    % make visible legend groups
    ghGray = hggroup; ghRed = hggroup;
    set(get(get(ghGray , 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    set(get(get(ghRed, 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    for ii = 1:numel(bins) % draw histogram bin by bin
        % Gray bin
        xl = bins(ii) - bw(ii)/2; xr = bins(ii);
        yd = 0; yu = hGray(ii);
        patch([xl xl; xl xr; xr xr], [yd yu; yu yu; yd yd], [0.5 0.5 0.5],...
            'EdgeColor','none', 'Parent', ghGray);
        hold on;
        % Red bins
        xl = bins(ii); xr = bins(ii) + bw(ii+1)/2;
        yd = 0; yu = hRed(ii);
        patch([xl xl; xl xr; xr xr], [yd yu; yu yu; yd yd], [1 0 0],...
            'EdgeColor','none', 'Parent', ghRed);
    end
end
if ~holdState
    hold off
end
end