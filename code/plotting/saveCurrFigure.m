function fname = saveCurrFigure(fname)
    ext = 'tiff';

    defaultDir = 'figures\';
    if isempty(strfind(fname, filesep)) && isempty(strfind(fname, '/'))
      fname = [defaultDir fname '.' ext];
    else
      fprintf('Saving figure to (%s)...\n', fname);
    end
      
    pathstr = fileparts(fname);
    if ~exist(pathstr, 'dir')
        fprintf('Creating new directory (%s)\n', pathstr);
        mkdir(pathstr);
    end
    %saveas(gcf, fname);
    %print(gcf, '-dtiff', fname);
    %    print(gcf, '-depsc', fname);
    % UNIX nodisplay mode
    switch ext
        case 'eps'
            driver_opts = {'-depsc','-painters'};
        case 'pdf'
            driver_opts = {'-dpdf','-painters','-r300'};
        otherwise
            driver_opts = {'-dtiff'};
    end

    % rand's secret sauce - maximize the screen figure
    scrsz = get(0,'ScreenSize');
    set(gcf, 'Position', [1 1 scrsz(3) scrsz(4)]);
    set(gcf, 'PaperPositionMode', 'auto');

    print(gcf, driver_opts{:}, fname);

    %fprintf('Converting to Pdf...\n');
    %sycmd = sprintf('mogrify -format pdf %s', fname);
    %system(sycmd);
end

