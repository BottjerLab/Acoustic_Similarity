function h1 = imagehist(A,varargin)
%imagehist(A,dotext, errors, subgrid)
%   Plot bivariate data as a color image with histograms along the margins.
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
[dotext, errors, subgrid] = optargs{:};

%% define grid on graph for subplots
Nfr = 2;
Nfc = 3;

xgrid = [0.025 0.175 0.3 0.85 0.875 0.925];
ygrid = [0.05 0.25 0.35 0.9];
xgrid = subgrid(1) + xgrid * subgrid(3);
ygrid = subgrid(2) + ygrid * subgrid(4);

for i = 1:Nfr
    for j = 1:Nfc
        p{Nfr+1-i,j}=[xgrid(2*j-1) ygrid(2*i-1) ...
            xgrid(2*j)-xgrid(2*j-1) ygrid(2*i)-ygrid(2*i-1)];
    end
end

%% get marginal distributions
nx = sum(A,2);
ny = sum(A,1);
cx = 1:size(A,2);
cy = 1:size(A,1);
%% get the common limits of the plots
xlims = [0 cx(end)]+0.5;
ylims = [0 cy(end)]+0.5;


%% plot marginal histogram along x direction
subplot('Position',p{2,2});

bar(cx,nx,1);
xlim(xlims);
set(gca,'XTickLabel','','YDir','reverse');
ylabel('sum');

%% plot marginal histogram along y direction
subplot('Position',p{1,1});
yoff = 0;
barh(cy-yoff,ny,1);
ylim(ylims-yoff);
set(gca,'YTickLabel', '', 'YDir','reverse', ...
    'XDir','reverse','YAxisLocation','right');
xlabel('sum');

%% plot color 2D
subplot('Position',p{1,2});
h1 = gca;
if min(min(A)) < 0
range = [min(min(A)) max(max(A))];
else
range = [0 max(max(A))];
end    
imagesc(A,range);
axis([xlims ylims]);
set(gca,'XTick',cx,'YTick',cy);

%% write text in color 2D image
if dotext
    for ii = 1:numel(cx)
        for jj = 1:numel(cy)
            
            if A(ii,jj)==ceil(A(ii,jj))
                s_mainnum = int2str(A(ii,jj));
            else
                s_mainnum = num2str(A(ii,jj),'%0.2f'); % two decimal places
            end
            if errors(ii,jj)>0
                str = sprintf('%s%s%0.2g',s_mainnum,setstr(177),errors(ii,jj));
            else
                str = s_mainnum;
            end
            hold on;
            h = text(jj, ii, str);
            hold off;
            set(h,'HorizontalAlignment','center');
            set(h,'VerticalAlignment','middle');
            set(h,'FontSize',12); % size should depend on grid size
        end
    end
end
% make colorbar
hcol=colorbar;
set(hcol,'Position',p{1,3});
