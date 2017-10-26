function plotHorzErrorBar(x, y, xwidth, col)
yh = 0.02*diff(ylim);
xx = [x-xwidth x+xwidth NaN x-xwidth x-xwidth NaN x+xwidth x+xwidth];
yy = [y y NaN y-yh y+yh NaN y-yh y+yh];
plot(xx,yy, '-','Color', col, 'LineWidth', 1.5);
%hold on;
%plot(x,y,'.','MarkerSize', 16, 'Color', col);
end