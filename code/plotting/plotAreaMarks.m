function hp = plotAreaMarks(regions, col, textLabel, yyOpt)

% regions is a structure
% col is a 3x1 RGB or 4x1 RGB + alpha triple
% textLabel is a bool indicating whether or not to display labels
% yyOpt is an optional variable designating 
if isempty(regions), 
    hp = NaN;
    return; 
end;
if nargin < 3 || ~islogical(textLabel)
    textLabel = false;
end
if nargin < 4 || isempty(yyOpt)
    yyOpt = ylim;
end

if nargin < 2 || isempty(col)
    col = [0.5 0.5 0.5]; %grey and opaque
end

% set alpha value to opaque by default
alphaVal = 1.0;
if numel(col) == 4
    alphaVal = col(4);
    col(4) = [];
end

nRegions = numel(regions);
starts = [regions.start]; stops = [regions.stop];
xx = [starts; starts; stops; stops];

yy = [yyOpt([1,2,2,1])]'; yy = yy(:,ones(1,nRegions));

if size(col,1) == 1 % single row, just one color provided
    vertexCols = col(ones(4*nRegions,1),:);
else % one color provided per patch
    assert(size(col,1) == nRegions);
    
    % an efficient way to copy each row 4 times
    % generate 1 1 1 1 2 2 2 2 3 3 3 3 etc....
    idxLookup = 1:nRegions; idxLookup = idxLookup(ones(4,1),:); idxLookup = idxLookup(:);
    vertexCols = col(idxLookup,:);
end
 
faceIdxs = reshape(1:4*nRegions,4,nRegions)';

% create the patch
hp = patch('Vertices',[xx(:) yy(:)], 'Faces', faceIdxs, ...
      'FaceColor','flat',...
      'FaceVertexCData', vertexCols,...
      'EdgeColor','flat','LineWidth',1,...
      'FaceAlpha',alphaVal);

 
% move lines in front of patches (FIXME)
cc = get(gca,'Children');
hdlIsLine = strcmp(get(cc,'Type'),'line');
set(gca,'Children', [cc(hdlIsLine); cc(~hdlIsLine)]);

if textLabel
     for ii = 1:nRegions
        thisLabel = regions(ii).type;
        if ~ischar(thisLabel), thisLabel = num2str(thisLabel); end
        xpos = (starts(ii) + stops(ii))/2;
        createTextLabel(xpos, thisLabel);
    end
end

end