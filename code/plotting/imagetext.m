function h1 = imagetext(A,varargin)
%imagehist(A,dotext, errors, subgrid)
%   Plot bivariate data as a color image
%   Returns a handle to the center plot
%   dotext = true, puts numbers in the plot
%
%   errors writes out errors with text
%
%   subgrid allows subplots within subplotting
%    where subgrid = [left bottom width height]
%   Example:
%
%      [x,y] = meshgrid(-10:10,-10:10);
%      gauss = exp(-(x.^2+y.^2)/25);
%      imagehist(gauss, false)

% note: doesn't do well with NaNs

% code from MATLAB blogs
numvarargs = length(varargin);
if numvarargs > 4
    error('imagehist:TooManyInputs', ...
        'requires at most 4 optional inputs');
end

% set defaults for optional inputs
optargs = {true, zeros(size(A)), [0 0 1 1], true};

% skip any new inputs if they are empty
newVals = ~cellfun('isempty', varargin);
% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);

% Place optional args in memorable variable names
[dotext, errors, subgrid, doColorbar] = optargs{:};

%% define grid on graph for subplots
xgrid = [0.13 0.885 0.91 0.955];
ygrid = [0.1 0.92];
%xgrid = [0.15 0.85]; 
Nfr = numel(ygrid)/2;
Nfc = numel(xgrid)/2;

xgrid = subgrid(1) + xgrid * subgrid(3);
ygrid = subgrid(2) + ygrid * subgrid(4);

for i = 1:Nfr
   for j = 1:Nfc
        p{Nfr+1-i,j}=[xgrid(2*j-1) ygrid(2*i-1) ...
                xgrid(2*j)-xgrid(2*j-1) ygrid(2*i)-ygrid(2*i-1)];
   end
end

%% get the common limits of the plots 
cx = 1:size(A,2);
cy = 1:size(A,1);
xlims = [0 cx(end)]+0.5;
ylims = [0 cy(end)]+0.5;

%% plot color 2D
subplot('Position',p{1,1});
h1 = gca;
if min(min(A)) < 0
    range = [min(A(:)) max(A(:))];
else
    range = [0 max(A(:))]; 
end

ih = imagesc(A,range);
%ih = imagesc(A,'AlphaData',~isnan(A),range); 
axis([xlims ylims]);
set(gca,'XTick',cx,'YTick',cy,'Color','none','Box','off','TickLength',[0 0]);

%% write text in color 2D image
if dotext
    for ii = 1:numel(cx)
        for jj = 1:numel(cy)
            
            entry = A(jj,ii);
            if isnan(entry) %don't print NaNs
                continue;
            end
            
            if entry==ceil(entry)
                s_mainnum = int2str(entry);
            else
                s_mainnum = num2str(entry,'%0.3f'); % two decimal places
            end
            
            if errors(jj,ii)>0
                % adds +- symbol, not very printable...
                str = sprintf('%s%s%0.2g',s_mainnum,setstr(177),errors(ii,jj)); 
            else
                str = s_mainnum;
            end
            
            hold on;
            h = text(ii,jj, str);
            hold off;
            %set(h, 'Color', [1 1 1]);
            set(h,'HorizontalAlignment','center');  
            set(h,'VerticalAlignment','middle'); 
            set(h,'FontSize',12); % size should depend on grid size
        end
    end
end

% make colorbar
if doColorbar
    hcol=colorbar;
    set(hcol,'Position',p{1,2});
end
