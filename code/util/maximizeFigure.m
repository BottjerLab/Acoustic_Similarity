function maximizeFigure
scrSize=get(0,'ScreenSize');
set(gcf,'Position', [1 1 scrSize(3:4)]); % assumes units are in pixels, which is default