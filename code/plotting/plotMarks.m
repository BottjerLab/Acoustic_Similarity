function plotMarks(regions, yc)
if isempty(regions), return; end;
if nargin < 2
    yc = ylim * [0; 1];
end
holdState = ishold;
hold on;

xac = [regions.start];
xbc = [regions.stop];
plot(xac,yc, 'cv',...
    'MarkerFaceColor','c', 'MarkerSize', 5);
plot(xbc,yc, 'mv',...
    'MarkerFaceColor','m', 'MarkerSize', 5);

xac = [xac; xac; NaN(size(xac))]; 
xac = xac(:);

xbc = [xbc; xbc; NaN(size(xbc))];
xbc = xbc(:);

yc = [ylim NaN]';
yc = yc(:,ones(numel(regions),1));
yc = yc(:);

if strcmp(get(gca,'YScale'),'log')
    semilogy(xac, yc,'c-',xbc, yc,'m-');
else
    plot(xac, yc,'c-',xbc, yc,'m-');
end

if ~holdState, hold off; end
end