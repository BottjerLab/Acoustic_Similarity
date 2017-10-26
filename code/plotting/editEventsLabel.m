function evsNew = editEventsLabel(evs,fs,doLabel)
% EDITEVENTSLABEL(EVS)
%   This function allows the user to edit event boundaries and labels.
%
%   The function plots a set of gray patches over each event in the active
%   figure.  Usually a waveform or some other line plot pertaining to
%   the signal should be behind it.  
%   The user can left-click and drag on events to adjust either the onset,
%   offset, or both.  The user can also right-click on events to relabel
%   them.  Clicking and dragging on a space without an event will create
%   a new event, while dragging the boundaries of an event past each other
%   results in deletion.
%
%   For the purpose of this 
%
% Known issues: 
% 1) Cursor may jump to the wrong side if patches are too close together
% 2) Can be slow if surface data is being displayed in same figure
% 3) If slow and you try to type a label too soon, focus moves to command window
% 4) Events cnnot be added when no events exist?
% Prereqs: active figure/axes that are appropriate to have marks
% evs is properly sorted

% TODO: link callbacks so that drags can be performed everywhere
% figure out a way to draw patches smartly over the spectrograms
% NB: for resizing to work properly, the axes property 'Units' should be
% normalized

% housekeeping, removing warning
RGBWarnID  = 'MATLAB:hg:patch:RGBColorDataNotSupported';
warnState = warning('query',RGBWarnID);
warning('off',RGBWarnID');
if nargin < 3
    doLabel = false;
end

% setting some defaults for sampling rate
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
newTitle = [oldTitle ' - Double click outside figure to exit and save, double click to play sound'];
if doLabel, newTitle = [newTitle ', right click to relabel']; end;
title(newTitle);

% hide any other patch handles that are being used here
otherPatchesOnAxis = findobj('Type','patch','Parent',gca);
set(otherPatchesOnAxis,'visible','off');

% create patch handle
if isempty(evs)
    % create a fake event and then delete it later
    hp = plotAreaMarks(initEvents(1),greycol);   
else
    hp = plotAreaMarks(evs,greycol);
end

% prepare handle for axes
hax = gca;

% give the label information to the patch handles
set(hp,'UserData',struct('labels', [evs.type]));
% transform labels into text objects
if doLabel
    drawLabelsInit(hp);
end

% some notes on action:
% when in drag mode, data gets added to the patch's userData field             
% when changing labels, the data gets changed in the patch's userData.labels field
set(gcf,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
set(gcf,'WindowButtonUpFcn',{@buttonUpFcn, [hax hp]});
set(gcf,'WindowButtonDownFcn',{@buttonDownFcn, [hax hp]});

% exit is triggered when mouseUp occurs outside the axis window
disp('Click outside to finish...');
waitfor(gcf,'WindowButtonMotionFcn','');
disp('Wrapping up...');
% clean up - restore any changed properties
title([oldTitle ' - Finishing']);
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

if doLabel
    % retrieve the labels
    labelHandles = getfield(get(hp,'UserData'),'labels');
    textLabels = get(labelHandles,'String');
    if ~iscell(textLabels), textLabels = {textLabels}; end;
    [evsNew.type] = textLabels{:};
    
    % get rid of the old text labels
    delete(labelHandles);
    set(otherPatchesOnAxis,'XData',get(hp,'XData'),'visible','on');
end
    
% get rid of the old patch
delete(hp);
% if we had an empty event structure to begin with, remove the placeholder
% first structure
if isempty(evs)
    evsNew(1) = [];
end
title(oldTitle);

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
    vertexColors(patchHover * 4 - 3,:) = hiCol;
    if ~isempty(lineHover)
        %vertexColors(patchHover * 4 + [-2,0],:) = lineCol(ones(2,1),:);
        currXData = get(hp,'XData');
        cursor(currXData(1+lineHover,patchHover),'on');
    else
        cursor(0,'off');
    end
else
    cursor(0,'off');
end

set(hp,'FaceVertexCData',vertexColors);
end

function cursor(xpos, status)
    lineHandle = findobj('Tag','cursor');
    if isempty(lineHandle)
        line([xpos xpos], ylim,'Color',lineCol,'Tag','cursor','LineWidth',1.5);
    elseif strcmp(status, 'off') && strcmp(get(lineHandle,'visible'), 'on')
        set(lineHandle,'visible','off');
    elseif strcmp(status, 'on')
        set(lineHandle,'XData',[xpos xpos],'visible',status);
    end
end

function buttonUpFcn(gcbo, eventdata, handles)
hax = handles(1); hp = handles(2);
% get point relative to window to determine if click lies outside
currPt = get(gcbo,'CurrentPoint'); currPt(1,1:2); 

oldUnits = get(gca,'Units');
set(gca,'Units','pixels');
axisWindow = get(hax,'Position');
set(gca,'Units',oldUnits);

% return cursor to original look
setptr(gcf,'arrow');

% exiting function - did we double-click outside the figure and 
% not as part of a drag?
if ~inRect(axisWindow, currPt) && ~isfield(get(hp,'UserData'),'lineHeld') && ...
    strcmp(get(gcbo,'SelectionType'),'open')

    % clearing the callbacks is the signal for the program to exit
    set(gcbo,'WindowButtonMotionFcn','');
    set(gcbo,'WindowButtonUpFcn','');
    set(gcbo,'WindowButtonDownFcn','');

elseif isfield(get(hp,'UserData'),'patchHeld') % finished clicking on an event
    % remove the data not pertaining to labels
    userData = get(hp,'UserData');
    set(hp,'UserData',struct('labels',userData.labels,'lastClicked',userData.patchHeld));
    
    %(1) negative intervals - delete
    %(2) overlapping intervals - merge    
    resolveOverlaps(hp);

    set(gcbo,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
else
    % clear last Clicked field
    userData = get(hp,'UserData');
    set(hp,'UserData',struct('labels',userData.labels)); % leave only the labels
end
    % links patches from other plots - a bit hacky, assumes all plots 
    % hold same patch pattern
    % note that  a naive linkprop doesn't work because we need the size of 
    % yData/colorData to vary dynamically
    % TODO: dynamically activate/deactivate linkprop
         otherPatches = findobj(gcbo, 'Type', 'patch');
         for ii = 1:numel(otherPatches)
             otherYData = get(otherPatches(ii),'YData');
             otherYData = otherYData(:,ones(1,size(get(hp,'YData'),2)));
             set(otherPatches(ii),'XData',get(hp,'XData'),...
                 'YData',otherYData,...
                 'FaceVertexCData',get(hp,'FaceVertexCData'));
         end
    
    % refresh labels
    if doLabel
        repositionLabels(hp);
    end
end

function buttonDownFcn(gcbo, eventdata, handles)
    hax = handles(1); hp = handles(2);
    currPt = get(hax,'CurrentPoint'); currPt = currPt(1,1:2);

    [patchClicked, lineClicked, hitWindow] = clickStatus(currPt, handles);
    if ~hitWindow, return; end;
    if doLabel && strcmp(get(gcbo,'SelectionType'),'alt') % if a right click, rename
        if ~isempty(patchClicked)
            textHandles = getfield(get(hp,'UserData'),'labels');
            set(textHandles(patchClicked),'Editing','on');
        end
    elseif strcmp(get(gcbo,'SelectionType'),'extend') % middle click
        disp('Debug Mode:');
        keyboard;
%         profile viewer;
%         profile on;
%         pause;
    elseif strcmp(get(gcbo,'SelectionType'),'open') % double click detection? play a sound
        % let people know we're working while we load the clip
%         set(gcf,'Pointer','watch');
%         disp('Clock on');
%         
        % get the clip data
        hline = findobj(hax,'Type','line');
        hline = hline(end); % it should be the one furthest back
        xWave = get(hline,'XData'); yWave = get(hline,'YData');
        if isempty(patchClicked) % what should we do if a patch is not clicked?
            % play the whole clip
            playSound(yWave, fs, true);
        else            
            % get the borders
            currXData = get(hp, 'XData');
            patchBorders = currXData(2:3,patchClicked);
            
            % cut the clip at the right points
            clipStart = find(xWave >= patchBorders(1),1);
            clipEnd = find(xWave >= patchBorders(2),1);
            clip = yWave(clipStart:clipEnd);
            
            % play the sound clip, while blocking
            playSound(clip, fs, true);
        end
        
        % ok, waiting's up
%         disp('Clock off');
%         set(gcf,'Pointer', 'arrow');
    elseif isempty(patchClicked) % clicked on an empty space
        % create new event
        currXData = get(hp, 'XData');
        currYData = get(hp, 'YData');
        currVertexColors = get(hp, 'FaceVertexCData');
        userData = get(hp,'UserData');
        % find where to insert new event
        if ~all(isnan(currXData))
            insertPt = find(currPt(1) <= [currXData(1,:) Inf], 1);
            currXData = [currXData(:,1:insertPt-1) currPt(1)*ones(4,1) currXData(:,insertPt:end)];
            currYData = currYData(:,[1 1:end]); %all columns are the same
            currVertexColors = currVertexColors([ones(1,4) 1:end], :);% we need four more rows, but the colors are all the same
            if doLabel
                newTextHandle = createTextLabel(currPt(1), ' ');
                userData.labels = [userData.labels(:,1:insertPt-1) newTextHandle userData.labels(:,insertPt:end)];
            end
        else
            % currXData is filled with a NaN box so that insertPt is not 
            % consistent with adding into an empty array 
            insertPt = 1;
            currXData = currPt(1) * ones(4,1);
            currYData = [ylim fliplr(ylim)]';
            if doLabel
                newTextHandle = createTextLabel(currPt(1), ' ');
                userData.labels=newTextHandle;  
            end
        end
        set(hp,'XData',currXData,'YData',currYData,'FaceVertexCData',currVertexColors);

        % handle dragging
        dragData = struct('lineHeld', [], ...
                          'startPt', currPt, ...
                          'patchHeld', insertPt, ...
                          'justCreated', true, ...
                          'origBounds', currPt(1)*ones(4,1), ...
                          'labels', userData.labels);
        
        % dragging data gets added to the patch userData              
        set(hp,'UserData',dragData);
        set(gcf,'WindowButtonMotionFcn',{@draggingFcn, [hax hp]});
        
        % set the cursor look
       setptr(gcf, 'fullcrosshair');
    else % we clicked on a patch
        xBounds = get(hp,'XData');
        userData = get(hp,'UserData');        
        dragData = struct('lineHeld', lineClicked, ...
                          'startPt', currPt, ...
                          'patchHeld', patchClicked,...
                          'justCreated', false,...
                          'origBounds', xBounds(:,patchClicked),...
                          'labels', userData.labels);
        
        set(hp,'UserData',dragData);
        set(gcf,'WindowButtonMotionFcn',{@draggingFcn, [hax hp]});
 
        % set cursor look
       if ~isempty(lineClicked)
           setptr(gcf,'fullcrosshair');
       else
           setptr(gcf,'hand');
       end
    end    
end

function draggingFcn(gcbo, eventdata, handles)
    hax = handles(1); hp = handles(2);
    currPt = get(hax,'CurrentPoint');
    
    currXData = get(hp,'XData');
    userData = get(hp, 'UserData');
    %nFaces = size(currXData,2);
    
    if ~isfield(userData,'patchHeld') || isempty(userData.patchHeld)
        error('editEvents:PatchNotClicked','Patch not Clicked, drag callback should not be set'); 
    end
        
    if userData.justCreated
        % if we just created an event, detect the drag motion
        if currPt(1) ~= userData.startPt(1)
            userData.justCreated = false;
            userData.lineHeld = 1 + (currPt(1) > userData.startPt(1));
        end
    end
    if ~isempty(userData.lineHeld) % moving one edge of the eventdata
        rIdxs = 2 * userData.lineHeld + [-1 0];
        currXData(rIdxs,userData.patchHeld) = currPt(1);
        
        set(hp,'XData',currXData);    
        cursor(currPt(1),'on');    

        % detect collisions immediately and quit drag
        if detectCollision(hp,userData.patchHeld,userData.lineHeld),
            set(gcbo,'WindowButtonMotionFcn',{@mouseOverFcn, [hax hp]});
            resolveOverlaps(hp);
            cursor(0,'off');
        end
    else % moving whole event
        currXData(:,userData.patchHeld) = userData.origBounds + currPt(1) - userData.startPt(1);        
        set(hp,'XData',currXData);    
    end
    
    if doLabel
        repositionLabels(hp);
    end
end

%%%%%%%% end callbacks %%%%%%%%%%%%%

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

function repositionLabels(hp)
    userData = get(hp,'UserData');
    currXData = get(hp,'XData');
    nFaces = size(currXData,2);
    for ii = 1:nFaces
        labelPos = get(userData.labels(ii),'Position');
        labelPos(1) = mean(currXData(:,ii));
        set(userData.labels(ii),'Position', labelPos);
    end
end

function resolveOverlaps(patchHandle)
    % cleans up intervals
    
    %keeping colors the same
    greyCol = [0.75 0.75 0.75];
    
    currXData = get(patchHandle,'XData');
    currYData = get(patchHandle,'YData');
    %currColData = get(patchHandle,'FaceVertexCData');
    
    if isempty(currXData), return; end; % nothing to resolve
    borders = currXData(2:3,:);
    userData = get(patchHandle,'UserData');
    
    % remove any negative-length intervals
    isNonPosLength = (borders(1,:) >= borders(2,:));
    borders(:,isNonPosLength) = [];
    
    % remove their label
    if doLabel
        delete(userData.labels(isNonPosLength));
        userData.labels(isNonPosLength) = [];
    end
    
    % merge overlapping regions 
    % since the number of regions is probably small (<10), we'll do this in a
    % naive way (better is with interval trees)
    ii = 1;
    nFaces = size(borders,2);
    while ii <= nFaces && nFaces > 1

        toMerge = find(borders(1,ii) >= borders(1,:) & borders(1,ii) <= borders(2,:) | ...
                       borders(2,ii) >= borders(1,:) & borders(2,ii) <= borders(2,:));
        
        if numel(toMerge) > 1
            borders(1,ii) = min(borders(1,toMerge));
            borders(2,ii) = max(borders(2,toMerge));
            toDelete = toMerge(toMerge ~= ii);
            
            % remove patch information
            borders(:,toDelete) = [];
            % remove labels
            if doLabel
                delete(userData.labels(toDelete));
                userData.labels(toDelete) = [];
            end
            % get resized # of patches
            nFaces = size(borders,2);
        else
            ii = ii + 1;
        end
    end
    
    if ~isempty(borders)
        currXData = borders([1 1 2 2],:);
        currYData = currYData(:,ones(1,nFaces)); % just copy the first row
        currColData = greyCol(ones(4*nFaces,1),:); % just copy the first color
    else
        currXData = NaN(4,1);
        currColData = greyCol(ones(4,1),:);
        currYData = NaN(4,1);
    end
   
    set(patchHandle,'XData',currXData,'YData',currYData,'FaceVertexCData',currColData,'UserData',userData);
    
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
    edgeFuzzFrac = 3e-3;
    edgeFuzz = diff(get(hax,'XLim')) * edgeFuzzFrac;

    yy = get(hax,'YLim');
    
    if currPt(2) > yy(2) || currPt(2) < yy(1), return; end;

    xBounds = get(hp,'XData'); 
    if isempty(xBounds), return; end; % nothing to click
    xBounds = xBounds(2:3,:);
    
    patchSeld = ...
        find(xBounds(1,:) - edgeFuzz <= currPt(1) & ...
        xBounds(2,:) + edgeFuzz >= currPt(1));

    if numel(patchSeld) > 1
        % find the one which is closer
        distsToCursor = min(xBounds(1,patchSeld) - currPt(1));
        [~,closest] = min(distsToCursor);
        patchSeld = patchSeld(closest);
    end
    if ~isempty(patchSeld)
        if     abs(xBounds(1,patchSeld) - currPt(1)) <= edgeFuzz, lineSeld = 1;
        elseif abs(xBounds(2,patchSeld) - currPt(1)) <= edgeFuzz, lineSeld = 2;
        end
    end
end

function drawLabelsInit(patchHandle)
    %prereq: userdata is already set with labels
    %converts labels to a set of txthandles
    currXData = get(patchHandle,'XData');
    borders = currXData(2:3,:);
    nFaces = size(borders,2);
    userData = get(patchHandle,'UserData');
    isString = ischar(userData.labels);
    txthandles = zeros(1,nFaces);
    for ii = 1:nFaces
        % convert to string if necessary
        thisLabel = userData.labels(ii);
        if ~isString, thisLabel = num2str(thisLabel); end
        txthandles(ii) = createTextLabel(mean(borders(:,ii)),thisLabel);
    end
    % labels is cell
    userData.labels = txthandles;
    set(patchHandle,'UserData',userData);
end
end
