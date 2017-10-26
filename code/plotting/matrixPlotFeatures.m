function matrixPlotFeatures(features, marks)
if nargin < 2
    marks = false(1,numel(features));
end
fields = fieldnames(features);
% remove all fields that are about variance
mean_fields = fields(cellfun(@isempty, strfind(fields,'_mean')));

N = numel(mean_fields);
logged = false(1,N);
for ii = 1:N
    % check for log scaling
    logged(ii) = any([features.(mean_fields{ii})] < 0);
end

nP = 5; % n per figure
for ii = 1:N
    for jj = 1:N
        figIdx = ceil(jj/nP) + (ceil(ii/nP) - 1) * ceil(N/nP);
        
        subII = mod(ii-1,nP)+1;
        subJJ = mod(jj-1,nP)+1;
        
        subPlotIdx = subJJ + (subII - 1) * nP; 
        
        figure(figIdx)
        hh = subplot(nP,nP, subPlotIdx);
        
        if ii == jj
            text(0,0.5,strrep(mean_fields{ii},'_','\_')); 
        else
            xx = [features.(mean_fields{jj})];
            yy = [features.(mean_fields{ii})];
            
            plot(xx(~marks), yy(~marks), 'r.',  ...
                 xx( marks), yy( marks), 'g.','MarkerSize',4);
            if ~logged(jj) 
                set(gca,'XScale','log');
            end
            if ~logged(ii) 
                set(gca,'YScale','log');
            end
            
        end

    end
end

end
