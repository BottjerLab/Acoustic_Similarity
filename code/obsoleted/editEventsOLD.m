function evsNew = editEvents(evs,fs)
% Prereqs: active figure/axes that are appropriate to have marks
% evs is properly sorted

% NB: for resizing to work properly, the axes property 'Units' should be
% normalized

% housekeeping, removing warning
RGBWarnID  = 'MATLAB:hg:patch:RGBColorDataNotSupported';
warnState = warning('query',RGBWarnID);
warning('off',RGBWarnID');

if isempty(evs), 
    evs = initEvents;
    if nargin < 2
        fs = 44100; % FIXME: a pure guess
        warning('editEvents:InputUninitialized','Events uninitialized, sampling rate may be incorrect...');
    end
else
    if nargin < 2
        fs = evs(1).idxStart/evs(1).start;
    end
end
greycol = [0.75 0.75 0.75];

% inform user of termination behavior
oldTitle = get(get(gca,'Title'),'String');
newTitle = 'Click outside figure to exit and save';
title(newTitle);

% create patch handle
if isempty(evs)
    % create a fake event and then delete it later
    hp = plotAreaMarks(initEvent,greycol);   
else
    hp = plotAreaMarks(evs,greycol);
end
% prepare handle for axes
hax = gca;

set(gcf,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
set(gcf,'WindowButtonUpFcn',{@buttonUpFcn, [hax hp]});
set(gcf,'WindowButtonDownFcn',{@buttonDownFcn, [hax hp]});

disp('Click outside to finish');

% exit is triggered when mouseUp occurs outside the axis window
waitfor(gcf,'WindowButtonMotionFcn','');

% clean up - restore any changed properties
title(oldTitle);
warning(warnState.state,RGBWarnID);

% read the events back from the edited patch handle

xdat = get(hp,'XData'); 
evsNew = initEvents(size(xdat,2));
if isempty(xdat), return; end;  % return an empty event structure if no marks

starts = num2cell(xdat(2,:)); idxStarts = num2cell(floor(xdat(2,:) * fs)); 
stops  = num2cell(xdat(3,:)); idxStops  = num2cell( ceil(xdat(3,:) * fs)); 
[evsNew.start] = starts{:}; [evsNew.stop] = stops{:};
[evsNew.idxStart] = idxStarts{:}; [evsNew.idxStop] = idxStops{:}; 
[evsNew.type] = deal(NaN);

% get rid of the old patch
delete(hp);

% if we had an empty event structure to begin with, remove the placeholder
% first structure
if isempty(evs)
    evsNew(1) = [];
end
end

%%%%%%%% begin callbacks %%%%%%%%%%%%%

function mouseOverFcn(gcbo, eventdata, handles)
% some default colors
greyCol = [0.75 0.75 0.75];
hiCol = [0.5 0.5 0.5];
lineCol = [0.8 0 0];

% unpack handles
hax = handles(1); hp = handles(2);

currPt = get(hax,'CurrentPoint'); currPt = currPt(1,1:2);

nFaces = size(get(hp,'XData'),2);
vertexColors = greyCol(ones(4 * nFaces,1),:);

[patchHover, lineHover] = clickStatus(currPt, handles);
if ~isempty(patchHover)
    % highlight that patch
    vertexColors(patchHover * 4 - 3,:) = hiCol; % what if numel(patchHover) > 1
    if ~isempty(lineHover)
        vertexColors(patchHover * 4 + [-2,0],:) = lineCol(ones(2,1),:);
    end
end

set(hp,'FaceVertexCData',vertexColors);
end

function buttonUpFcn(gcbo, eventdata, handles)
hax = handles(1); hp = handles(2);
% get point relative to window to determine if click lies outside
currPt = get(gcbo,'CurrentPoint'); currPt(1,1:2); 

oldUnits = get(gca,'Units');
set(gca,'Units','pixels');
axisWindow = get(hax,'Position');
set(gca,'Units',oldUnits);

% exiting function - did we click outside the figure and not as part of a drag?
if ~inRect(axisWindow, currPt) && isempty(get(hp,'UserData'));

    % clearing the callbacks is the signal for the program to exit
    set(gcbo,'WindowButtonMotionFcn','');
    set(gcbo,'WindowButtonUpFcn','');
    set(gcbo,'WindowButtonDownFcn','');

elseif ~isempty(get(hp,'UserData')) % finished editing drag
    % clean up intervals

    %(1) negative intervals - delete
    %(2) overlapping intervals - merge    
    resolveOverlaps(hp);

    set(hp,'UserData','');    
    set(gcbo,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
    disp('Letting go');
else
    disp('Still holding on...');
end
end

function buttonDownFcn(gcbo, eventdata, handles)
    hax = handles(1); hp = handles(2);
    currPt = get(hax,'CurrentPoint'); currPt = currPt(1,1:2);

    [patchClicked, lineClicked, hitWindow] = clickStatus(currPt, handles);
    if ~hitWindow, return; end;
    if isempty(patchClicked) % nothing, or create new event?
        % create new event
        currXData = get(hp, 'XData');
        currYData = get(hp, 'YData');
        currVertexColors = get(hp, 'FaceVertexCData');
        % find where to insert new event
        insertPt = find(currPt(1) <= [currXData(1,:) Inf], 1);
        currXData = [currXData(:,1:insertPt-1) currPt(1)*ones(4,1) currXData(:,insertPt:end)];
        currYData = currYData(:,[1 1:end]); %all columns are the same
        currVertexColors = currVertexColors([ones(1,4) 1:end], :);% we need four more rows, but the colors are all the same
        
        disp('Creating patch');
        set(hp,'XData',currXData,'YData',currYData,'FaceVertexCData',currVertexColors);

        % handle dragging 
        dragData = struct('lineHeld', [], ...
                          'startPt', currPt, ...
                          'patchHeld', insertPt, ...
                          'justCreated', true, ...
                          'origBounds', []);
        set(hp,'UserData',dragData);
        set(gcf,'WindowButtonMotionFcn',{@draggingFcn, [hax hp]});
    else
        xBounds = get(hp,'XData');
        
        dragData = struct('lineHeld', lineClicked, ...
                          'startPt', currPt, ...
                          'patchHeld', patchClicked,...
                          'justCreated', false,...
                          'origBounds', xBounds(:,patchClicked));
                      
        
        set(hp,'UserData',dragData);
        set(gcf,'WindowButtonMotionFcn',{@draggingFcn, [hax hp]});
    end
end

function draggingFcn(gcbo, eventdata, handles)
    hax = handles(1); hp = handles(2);
    currPt = get(hax,'CurrentPoint');
    
    currXData = get(hp,'XData');
    userData = get(hp, 'UserData');
    %nFaces = size(currXData,2);
    
    if isempty(userData.patchHeld)
        error('editEvents:PatchNotClicked','Patch not Clicked, drag callback should not be set'); 
    end
        
    if userData.justCreated
        % if we just created an event, detect the drag motion
        if currPt(1) ~= userData.startPt(1)
            userData.justCreated = false;
            userData.lineHeld = 1 + (currPt(1) > userData.startPt(1));
        end
        set(hp,'UserData',userData);
    end
    if ~isempty(userData.lineHeld) % moving one edge of the eventdata
        rIdxs = 2 * userData.lineHeld + [-1 0];
        currXData(rIdxs,userData.patchHeld) = currPt(1);
        
        set(hp,'XData',currXData);    

        % detect collisions immediately and quit drag
        if detectCollision(hp,userData.patchHeld,userData.lineHeld),
            set(gcbo,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
            resolveOverlaps(hp);
        end
    else % moving whole event
        currXData(:,userData.patchHeld) = userData.origBounds + currPt(1) - userData.startPt(1);        
        set(hp,'XData',currXData);    
    end
    
end

%%%%%%%% end callbacks %%%%%%%%%%%%%

function resolveOverlaps(patchHandle)
    % cleaning up intervals, 
    
    %keeping colors the same
    greyCol = [0.75 0.75 0.75];
    
    currXData = get(patchHandle,'XData');
    currYData = get(patchHandle,'YData');
    %currColData = get(patchHandle,'FaceVertexCData');
    
    if isempty(currXData), return; end; % nothing to resolve
    borders = currXData(2:3,:);

    % remove any negative-length intervals
    isNonPosLength = (borders(1,:) >= borders(2,:));
    borders(:,isNonPosLength) = [];
    
    % merge overlapping regions 
    % since the number of regions is probably small (<10), we'll do this in a
    % naive way (better is with interval trees)
    ii = 1;
    nFaces = size(borders,2);
    while ii <= nFaces && nFaces > 1

        toMerge = find(borders(1,ii) >= borders(1,:) & borders(1,ii) <= borders(2,:) | ...
                       borders(2,ii) >= borders(1,:) & borders(2,ii) <= borders(2,:));
        
        if numel(toMerge) > 1
            disp('Merging...')
            borders(1,ii) = min(borders(1,toMerge));
            borders(2,ii) = max(borders(2,toMerge));
            toDelete = toMerge(toMerge ~= ii);
            borders(:,toDelete) = [];
            nFaces = size(borders,2);
        else
            ii = ii + 1;
        end
    end
    
    currXData = borders([1 1 2 2],:);
    currYData = currYData(:,ones(1,nFaces)); % just copy the first row
    currColData = greyCol(ones(4*nFaces,1),:); % just copy the first color
    set(patchHandle,'XData',currXData,'YData',currYData,'FaceVertexCData',currColData);
end

function didCollide = detectCollision(patchHandle,activePatch, boundarySide)
    currXData = get(patchHandle,'XData');
    if isempty(currXData), return; end; % nothing to collide 
    borders = currXData(2:3,:);
    nFaces = size(borders,2);
    
    adjBorder = NaN; 
    if activePatch > 1 && boundarySide == 1
        adjBorder = borders(2,activePatch - 1);
    elseif activePatch < nFaces && boundarySide == 2
        adjBorder = borders(1,activePatch + 1);
    end

    didCollide = (borders(1,activePatch) >= borders(2,activePatch)) || ...
        borders(boundarySide,activePatch) <= adjBorder && boundarySide == 1 || ...
        borders(boundarySide,activePatch) >= adjBorder && boundarySide == 2;
end

function foo = inRect(win, pt)
foo = win(1) <= pt(1) && pt(1) < win(1) + win(3) && ...
    win(2) <= pt(2) && pt(2) < win(2) + win(4);
end

function [patchSeld, lineSeld, hitWindow] = clickStatus(currPt, handles)
% returns empties on default
    patchSeld = []; lineSeld = [];
    hax = handles(1); hp = handles(2);
    
    win([1 3]) = get(hax,'XLim'); win(3) = win(3) - win(1);
    win([2 4]) = get(hax,'YLim'); win(4) = win(4) - win(2);
    hitWindow = inRect(win, currPt); if ~hitWindow, return, end;
    
    % how 'fat' should our edge be for us to highlight/grab it?
    edgeFuzzFrac = 0.0035;
    edgeFuzz = diff(get(hax,'XLim')) * edgeFuzzFrac;

    yy = get(hax,'YLim');
    
    if currPt(2) > yy(2) || currPt(2) < yy(1), return; end;

    xBounds = get(hp,'XData'); 
    if isempty(xBounds), return; end; % nothing to click
    xBounds = xBounds(2:3,:);
    
    patchSeld = ...
        find(xBounds(1,:) - edgeFuzz <= currPt(1) & ...
        xBounds(2,:) + edgeFuzz >= currPt(1),1);

    if     any(abs(xBounds(1,:) - currPt(1)) <= edgeFuzz), lineSeld = 1;
    elseif any(abs(xBounds(2,:) - currPt(1)) <= edgeFuzz), lineSeld = 2;
    end
end