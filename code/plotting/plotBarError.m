function handles_bars=plotBarError(vals,err, ylims, sym, colors)
% plot bar graph with errors
% vals and err are NbarsxNgroups, the values and SEMs, respectively
% sym is the statistics marker indicating significance
Nbars = size(vals,1);
Ngroups = size(vals,2);
pats = {'*','o','+'};
if nargin<2 || isempty(err)
    err=zeros(size(vals));
end
autoset = (nargin<3) || isempty(ylims);
if nargin<4 || isempty(sym) %result of t-tests
    sym=zeros(size(vals));
end
if nargin < 5 || isempty(colors)
    colors = jet(Ngroups);
end
nargin;

topFlag = 1;

handles_bars=bar(vals,'grouped');
set(gca,'XTick',1:Nbars);
hold on;

bar_ratio = 1;
bar_abscissa=zeros(Nbars,Ngroups);
bar_width = zeros(Nbars,Ngroups);
for ii = 1:Ngroups
    % dive into bar configs to get x coords
    bar_edges = get(get(handles_bars(ii),'Children'),'Xdata');
    bar_abscissa(:,ii) = mean(bar_edges([1 3],:))';
    bar_width(:,ii) = bar_edges(3,:)-bar_edges(1,:);
    
    % assign colors to each bar group here
    set(handles_bars(ii), 'FaceColor',colors(ii,:));
end

% draw errorbars
for ii = 1:Ngroups
    if ~topFlag
        errorbar(bar_abscissa(:,ii), vals(:,ii), err(:,ii), 'k', ...
                 'linestyle', 'none', 'linewidth', 2);
    else
        xcoord = [bar_abscissa(:,ii) bar_abscissa(:,ii) ...
                  bar_abscissa(:,ii) + bar_ratio/2*bar_width(:,ii) ...
                  bar_abscissa(:,ii) - bar_ratio/2*bar_width(:,ii) ...
                  NaN(Nbars,1)]';
              
        yErrSign = sign(vals(:,ii));
        yHBarPos = vals(:,ii) + yErrSign.* err(:,ii);
        ycoord = [vals(:,ii) yHBarPos yHBarPos yHBarPos ...
                  NaN(Nbars,1)]';
        xcoord = xcoord(:);
        ycoord = ycoord(:);
        plot(xcoord, ycoord,'k','HandleVisibility','off');
    end
    
end

if ~autoset
    ylim(ylims);
end
if any(sym>0)
    ymax = ylim;
    ymax = ymax(2);
    for isym = 1:max(sym)
        beats = find(sym==isym);
        maxy = max(vals(:),[],2);
        ydots = (maxy(beats)+ymax)/2;
        plot(beats,ydots,['k' pats{isym}]);
    end
end
hold off;
xlim([0 Nbars] + 0.5);
end
