function handle = createTextLabel(xpos, label)
% label must be a string
handle = text(xpos, ylim * [0.05 0.95]', label,...
                'Color','white','BackgroundColor','red',...
                'HorizontalAlignment','center','FontWeight','bold');
end