function fig = plotLineTransparentSEM(x, y, SEM, col, varargin)
% plotLineTransparentSEM
% plots x and y, 1D vectors,  with color specified by an RGB row vector
% convert to row vectors
if size(x,1) > 1 
    if size(x,2) > 1
        error('x must be a vector.');
    else
        x = x';
    end
end
if size(y,1) > 1 
    if size(y,2) > 1
        error('y must be a vector.');
    else
        y = y';
    end
end
if size(SEM,1) > 1 
    if size(SEM,2) > 1
        error('SEM must be a vector.');
    else
        SEM = SEM';
    end
end


if ~(all(size(x) == size(y)) && all(size(x) == size(SEM)))
    size(x)
    size(y)
    size(SEM)
    error('x, y, and SEM must be the same size.');
end

if nargin < 4
    error('Function requires a RGB row vector as a color.');
end

washeld = ishold;
if ~isempty(varargin)
    plot(x,y,'-','color',col);
    %    plot(x,y,'-','color',col,varargin{:});
else
    plot(x,y,'-','color',col);
end    

% fill lines
hold on;
if usejava('desktop')% && ~strcmp(version, '8.1.0.604 (R2013a)') % only desktop versions of matlab can print transparency
    hfill = fill([x fliplr(x)], [y+SEM fliplr(y-SEM)], col, ...
                 'edgecolor', 'none'); % plot region
    alpha(0.2); % set transparency
    set(gcf,'Renderer', 'OpenGL'); % necessary for transparency, but won't do log scales =p (do the logs/customize tick marks yourself if you want log scale)
    set(hfill,'HandleVisibility','off'); % remove legend entry
else
    lighten = 0.4;
    lightCol = col * (1 - lighten) + lighten;
    plot(x, y+SEM, '-', 'color', lightCol, 'HandleVisibility','off', varargin{:});
    plot(x, y-SEM, '-', 'color', lightCol, 'HandleVisibility','off', varargin{:});
end
if ~washeld, hold off; end
    
fig = gca;
end

