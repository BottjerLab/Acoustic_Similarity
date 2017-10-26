function [hcomponent, hcontainer] = javacomponentundoc_helper(varargin)
% Copyright 2010-2012 The MathWorks, Inc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old javacomponent implementation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin>=1)
    component = varargin{1};
end
if (nargin>=2)
    position = varargin{2};
end

if (nargin>=3)
    parent = varargin{3};
end

if (nargin==4)
    callback =  varargin{4};
end

if (usejavacomponent == 0)
    error(message('MATLAB:javacomponent:FeatureNotSupported'));
end

if ~isempty(nargchk(1,4,nargin))  %#ok
    error('MATLAB:javacomponent',getString(message('MATLAB:javacomponent:IncorrectUsage')));
end

if nargin < 4
    callback = '';
else
    if ~iscell(callback)
        error('MATLAB:javacomponent',getString(message('MATLAB:javacomponent:IncorrectUsage')));
    end
end

if nargin < 3
    parent = gcf;
end

if nargin < 2
    position = [20 20 60 20];
end

parentIsFigure = false;
hParent = handle(parent);
% g500548 - changed to use ishghandle.
if ( ishghandle(hParent, 'figure') || ...
        ishghandle(hParent, 'uicontainer') || ...
        ishghandle(hParent, 'uiflowcontainer') || ...
        ishghandle(hParent, 'uigridcontainer'))
    parentIsFigure = true;
    
    peer = getJavaFrame(ancestor(parent,'figure'));
elseif ishghandle(hParent, 'uitoolbar')
    peer = get(parent,'JavaContainer');
    if isempty(peer)
        drawnow;
        peer = get(parent,'JavaContainer');
    end
elseif (ishghandle(hParent, 'uisplittool') || ...
        ishghandle(hParent, 'uitogglesplittool'))
    parPeer = get(get(hParent,'Parent'),'JavaContainer');
    if isempty(parPeer)
        drawnow;
    end
    peer = get(parent,'JavaContainer');
else
    error(message('MATLAB:javacomponent:InvalidParentHandle', getString(message('MATLAB:javacomponent:IncorrectUsage'))))
end

if isempty(peer)
    error(message('MATLAB:javacomponent:JavaFigsNotEnabled'))
end

hUicontainer = [];
hgp = [];
returnContainer = true;

if ischar(component)
    % create from class name
    component = javaObjectEDT(component);
elseif iscell(component)
    % create from class name and input args
    component = javaObjectEDT(component{:});
elseif ~isa(component,'com.mathworks.hg.peer.FigureChild')
    % tag existing object for auto-delegation - unless it is a FigureChild
    javaObjectEDT(component);
end

% Promote the component to a handle object first. It seems once a java
% object is cast to a handle, you cannot get another handle with
% 'callbackproperties'.
if ~isjava(component)
    component = java(component);
end
hcomponent  = handle(component,'callbackProperties');

if nargin == 1
    hgp = handle(peer.addchild(component));
    % parent must be a figure, we default to gcf upstairs
    createPanel;
    hgp.setUIContainer(double(hUicontainer));
else
    if parentIsFigure
        if isnumeric(position)
            if isempty(position)
                position = [20 20 60 20];
            end
            % numeric position is not set here, rely on the uicontainer
            % listeners below.
            hgp = handle(peer.addchild(component));
            createPanel;
            hgp.setUIContainer(double(hUicontainer));
        elseif ...
                isequal(char(position),char(java.awt.BorderLayout.NORTH)) || ...
                isequal(char(position),char(java.awt.BorderLayout.SOUTH)) || ...
                isequal(char(position),char(java.awt.BorderLayout.EAST))  || ...
                isequal(char(position),char(java.awt.BorderLayout.WEST))  || ...
                isequal(char(position),char('Overlay'))
            hgp = handle(peer.addchild(component, position));
            returnContainer = false;
        else
            error(message('MATLAB:javacomponent:InvalidPosition', getString(message('MATLAB:javacomponent:IncorrectUsage'))))
        end
    else
        % Adding component to the toolbar.
        % component position is ignored
        peer.add(component);
        hUicontainer = parent; % toolbar.
        handles = getappdata(hUicontainer, 'childhandles');
        handles = [handles, hcomponent];
        setappdata(hUicontainer, 'childhandles', handles);
    end
    
    % make sure the component is on the screen so the
    % caller can interact with it right now.
    % drawnow;
end

if returnContainer
    configureComponent();
end

% If asked for callbacks, add them now.
if ~isempty(callback)
    % The hg panel is the best place to store the listeners so they get
    % cleaned up asap. We can't do that if the parent is a uitoolbar so we
    % just put them on the toolbar itself.
    lsnrParent = hgp;
    if isempty(lsnrParent)
        lsnrParent = hParent;
    end
    if mod(length(callback),2)
        error('MATLAB:javacomponent',usage);
    end
    for i = 1:2:length(callback)
        lsnrs = getappdata(lsnrParent,'JavaComponentListeners');
        l = javalistener(component, callback{i}, callback{i+1});
        setappdata(lsnrParent,'JavaComponentListeners',[l lsnrs]);
    end
end

if returnContainer
    hcontainer = hUicontainer;
else
    hcontainer = [];
end

    function createPanel
        % add delete listener
        hUicontainer = hgjavacomponent('Parent',parent,'Units', 'Pixels','Serializable','off');
        
        set(hUicontainer, 'UserData', char(component.getClass.getName)); % For findobj queries.
        if isa(java(hgp), 'com.mathworks.hg.peer.FigureChild')
            set(hUicontainer, 'FigureChild', hgp);
        end
        if isa(java(hcomponent), 'javax.swing.JComponent')
            % force component to be opaque if it's a JComponent. This prevents
            % lightweight component from showing the figure background (which
            % may never get a paint event)
            hcomponent.setOpaque(true);
        end
        set(hUicontainer, 'JavaPeer', hcomponent);
        
        if returnContainer
            % add move/resize listener to the hgjavacomponent
            addlistener(hUicontainer, 'PixelBounds', 'PostSet', @handleResize);
            
            % add visible listener
            addlistener(hUicontainer, 'Visible', 'PostSet', @handleVisible);
            
            %Parent was set before we get here. Hence we need to explicitly
            %walk up and attach listeners. For subsequent parent changes,
            %the parent property postset listener callback will take care of setting
            %up the hierarchy listeners to listen for position and visible changes
            createHierarchyListeners(hUicontainer, @handleVisible);
            
            % add parent listener
            addlistener(hUicontainer, 'Parent', 'PreSet', @handlePreParent);
            
            % add parent listener
            addlistener(hUicontainer, 'Parent', 'PostSet', @(o,e) handlePostParent(o,e,@handleVisible));
            
            % force though 1st resize event
            set(hUicontainer,'Position', position);
        else
            % For the BorderLayout components, we don't really want the
            % hUicontainer to show. But, it provides a nice place for cleanup.
            set(hUicontainer,'Visible', 'off', 'HandleVisibility', 'off');
            % Set position out of the figure to work around a current bug
            % due to which invisible uicontainers show up when renderer is
            % OpenGL (G.
            set(hUicontainer, 'Position', [1 1 0.01 0.01]);
        end
        
        if isa(component,'com.mathworks.hg.peer.FigureChild')
            component.setUIContainer(double(hUicontainer));
        end
        
        %Handles move or resize of the hgjavacomponent.
        function handleResize(obj, evd) %#ok - mlint
            %Deal with early warnings
            [lastWarnMsg, lastWarnId] = lastwarn;
            oldPBWarning = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:PixelBounds');
            pb = get(hUicontainer,'PixelBounds');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % HACK FOR HMI - When X-location is negative, we are
            % messing up the width of the component because we are
            % going from pixel bounds (xmin,ymin,xmax,ymax) to actual
            % bounds on the screen (width and height calculated). There is 
            % a floor/ceil issue with the X-min not quite adjusted by
            % the X-max so as to keep the width correct. I am forcing it
            % here. An ideal fix would have been in the pixelbounds
            % calculation but it seems like too much work in HG1 and hence
            % I am using a narrow scoped fix for now  - 815546
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            jcPos = getpixelposition(hUicontainer,true);
            if (jcPos(1) < 0)
                pb(3) = pb(1) + jcPos(3);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hgp.setPixelBounds(pb);
            warning(oldPBWarning.state, 'MATLAB:HandleGraphics:ObsoletedProperty:PixelBounds');
            % restore the last warning thrown
            lastwarn(lastWarnMsg, lastWarnId);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % g482174 - Bandage solution to support hierarchy visibility changes
        % We will look for any invisible parent container to see if we can
        % show the javacomponent or not.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function handleVisible(obj, evd) %#ok - mlint
            source = evd.AffectedObject;
            if ishghandle(source,'hgjavacomponent')
                hgp.setVisible(strcmp(get(source,'Visible'),'on'));
            else
                if (strcmp(get(source,'Visible'),'off'))
                    setInternalVisible(hUicontainer, component, false);
                else
                    setInternalVisible(hUicontainer, component, isVisibleInCurrentLineage(hUicontainer));
                end
            end
            drawnow update;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % isVisibleInCurrentLineage - returns whether the javacomponent can
        % be visible or not in its current hierarchy.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function visible = isVisibleInCurrentLineage(hUicontainer)
            stParent = get(hUicontainer,'Parent');
            hgjcompVisible = strcmp(get(hUicontainer,'Visible'),'on');
            visible = hgjcompVisible && strcmp(get(stParent,'Visible'),'on');
            while (~ishghandle(stParent,'root') && visible)
                stParent = get(stParent,'Parent');
                visible = visible && strcmp(get(stParent,'Visible'),'on');
            end
        end
        
        
        
        function handlePreParent(obj, evd) %#ok - mlint
            oldfig = ancestor(hUicontainer, 'figure');
            removecomponent = true;
            % NewValue field is absent in MCOS and hence we need to do the
            % following safely.
            if ~isempty(findprop(evd,'NewValue'))
                newfig = ancestor(evd.NewValue, 'figure');
                removecomponent = ~isempty(newfig) && ~isequal(oldfig, newfig);
            end
            %We are losing on this optimization(event may not have NewValue). We always
            %remove and add upon reparenting. We do not have the context of
            %the new parent in the preset to do a compare to see if we are
            %being parented to the same parent again. We hope that this is
            %not done often.
            if  (removecomponent)
                peer = getJavaFrame(oldfig);
                peer.remove(component);
            end
        end
        
        function handlePostParent(obj, evd, visibleCbk) %#ok - mlint
            createHierarchyListeners(hUicontainer, visibleCbk);
            oldfig = ancestor(parent, 'figure');
            newfig = ancestor(evd.AffectedObject,'figure');
            addcomponent = true;
            %Before, we could decide if we want to re-add the
            %javacomponent. Now we have to always add.
            if  ~isempty(findprop(evd,'NewValue'))
                newfig = ancestor(evd.NewValue, 'figure');
                addcomponent = ~isempty(newfig) && ~isequal(oldfig, newfig);
            end
            
            if addcomponent
                peer = getJavaFrame(newfig);
                hgp= handle(peer.addchild(component));
                %When we reparent ourself (javacomponent), the truth about
                %whether we are visible or not needs to be queried from
                %the proxy and its current lineage. Change visibility
                %without changing state on the hUicontainer
                setInternalVisible(hUicontainer, component, isVisibleInCurrentLineage(hUicontainer));
                if isa(java(hgp), 'com.mathworks.hg.peer.FigureChild')
                    % used by the uicontainer C-code
                    setappdata(hUicontainer, 'FigureChild', java(hgp));
                end
                parent = newfig;
            end
            handleResize([],[]);
        end
    end

    function configureComponent
        set(hUicontainer,'DeleteFcn', {@containerDelete, hcomponent});
        %         addlistener(java(hcomponent), 'ObjectBeingDestroyed', @(o,e)componentDelete(o,e,hUicontainer, parentIsFigure));
        temp = handle.listener(hcomponent, 'ObjectBeingDestroyed', @(o,e)componentDelete(o,e,hUicontainer, parentIsFigure));
        save__listener__(hcomponent,temp);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g637916 - Visibility of the hgjavacomponent is posing issues due to the
% fact that the state is being used to control visibility and hence when
% the parent's are turned visible off, we need an alternate api to make it
% go away from the screen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setInternalVisible(hUicontainer, component ,vis)
cleaner = warnSuppressAndSetInternalVisible(hUicontainer, component ,vis);
delete(cleaner);
end

function cleaner = warnSuppressAndSetInternalVisible(hUicontainer, component ,vis)
[ state.lastWarnMsg, state.lastWarnId ] = lastwarn;
state.usagewarning = warning('off','MATLAB:hg:javacomponent');
cleaner = onCleanup(@() warnRestore(state));
setVisibility(handle(hUicontainer),vis);
if (isa(component,'com.mathworks.hg.peer.FigureChild'))
    component = component.getFigureComponent;
end
setVisible(component, vis);
end

function warnRestore(state)
lastwarn(state.lastWarnMsg, state.lastWarnId);
warning(state.usagewarning);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function containerDelete(obj, evd, hc) %#ok - mlint
obj = handle(obj);
if ishghandle(handle(obj), 'uitoolbar') || ...
        ishghandle(handle(obj),'uisplittool') || ...
        ishghandle(handle(obj),'uitogglesplittool');
    childHandles = getappdata(obj, 'childhandles');
    delete(childHandles(ishandle(childHandles)));
else
    if ishandle(hc)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We are giving ourselves a hook to run resource cleanup functions
        % like freeing up callbacks. This is important for uitree because
        % the expand and selectionchange callbacks need to be freed when
        % the figure is destroyed. See G769077 for more information.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        removeJavaCallbacks(hc);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(hc);
    end
end
end


function componentDelete(obj, evd, hUicontainer, parentIsFigure) %#ok - mlint
if (parentIsFigure)
    % This java component is always deleted before hUicontainer. It is
    % ensured by calling component deletion in function containerDelete.
    % hUicontainer becomes invalid when delete(hUicontainer) below is run.
    parent = ancestor(hUicontainer,'figure');
    
    peer = getJavaFrame(parent);
    
    if any(ishandle(obj))
        removeobj = java(obj);
        if ~isempty(get(hUicontainer,'FigureChild'))
            removeobj  = get(hUicontainer,'FigureChild');
        end
        peer.remove(removeobj);
    end
    
    % delete container if it exists
    if any(ishghandle(hUicontainer))
        delete(hUicontainer);
    end
else
    parent = hUicontainer; % toolbar or split tool
    if ~ishandle(parent) || ~ishandle(obj)
        % The toolbar parent or the component has been deleted. Bail out.
        % Toolbar clears all javacomponents after itself.
        return;
    end
    
    % For uisplittool and uitogglesplittool objects
    % The parent may have done this deletion for us
    % already.
    hPar = get(parent,'Parent');
    if ishghandle(handle(hPar),'uitoolbar')
        parPeer = get(hPar,'JavaContainer');
        if isempty(parPeer)
            return;
        end
    end
    
    peer = get(parent, 'JavaContainer');
    if ~isempty(peer)
        peer.remove(java(obj));
    end
end
end

function hdl=javalistener(jobj, eventName, response)
try
    jobj = java(jobj);
catch ex  %#ok
end

% make sure we have a Java objects
if ~ishandle(jobj) || ~isjava(jobj)
    error(message('MATLAB:javacomponent:invalidinput'))
end

hSrc = handle(jobj,'callbackproperties');
allfields = sortrows(fields(set(hSrc)));
alltypes = cell(length(allfields),1);
j = 1;
for i = 1:length(allfields)
    fn = allfields{i};
    if ~isempty(strfind(fn,'Callback'))
        fn = strrep(fn,'Callback','');
        alltypes{j} = fn;
        j = j + 1;
    end
end
alltypes = alltypes(~cellfun('isempty',alltypes));

if nargin == 1
    % show or return the possible events
    if nargout
        hdl = alltypes;
    else
        disp(alltypes)
    end
    return;
end

% validate event name
valid = any(cellfun(@(x) isequal(x,eventName), alltypes));

if ~valid
    error(message('MATLAB:javacomponent:invalidevent', class( jobj ), char( cellfun( @(x) sprintf( '\t%s', x ), alltypes, 'UniformOutput', false ) )'))
end

hdl = handle.listener(handle(jobj), eventName, ...
    @(o,e) cbBridge(o,e,response));
    function cbBridge(o,e,response)
        hgfeval(response, java(o), e.JavaEvent)
    end
end

function createHierarchyListeners(hUicontainer, visCbk)
deleteExistingHierarchyListeners(hUicontainer, []);
hUicontainer = handle(hUicontainer);
parent = get(hUicontainer,'Parent');
% Walk up instance hierarchy and put a listener on all the
% containers. We don't need a listener on the figure.
while ~ishghandle(parent,'root')
    %Set up all the visible listeners
    createVisibleListener(parent, hUicontainer, visCbk);
    %Keep walking up
    parent = get(parent,'Parent');
end

% When the hgjavacomponent goes away, clean all listeners
addlistener(hUicontainer, 'ObjectBeingDestroyed', @(o,e) deleteExistingHierarchyListeners(o,e));
end

function deleteExistingHierarchyListeners(src,~)
hUicontainer = handle(src);

%Delete all the visibility listeners
if isListenerData(hUicontainer, 'VisiblilityListeners')
    appdata = getListenerData(hUicontainer, 'VisiblilityListeners');
    cellfun(@(x) delete(x), appdata, 'UniformOutput',false);
    setListenerData(hUicontainer, 'VisiblilityListeners',{});
end



end

function createVisibleListener(object, hUicontainer, visCbk)
%Attach visibility listeners
visContainerListnr = addlistener(object, 'Visible','PostSet', visCbk);
if (~isListenerData(hUicontainer, 'VisiblilityListeners'))
    setListenerData(hUicontainer, 'VisiblilityListeners',{});
end
visibilityAppdata = getListenerData(hUicontainer, 'VisiblilityListeners');
visibilityAppdata{end+1} = visContainerListnr;
setListenerData(hUicontainer, 'VisiblilityListeners', visibilityAppdata);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function javaFrame = getJavaFrame(f)
% store the last warning thrown
[ lastWarnMsg, lastWarnId ] = lastwarn;

% disable the warning when using the 'JavaFrame' property
% this is a temporary solution
oldJFWarning = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
javaFrame = get(f,'JavaFrame');
warning(oldJFWarning.state, 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

% restore the last warning thrown
lastwarn(lastWarnMsg, lastWarnId);
end
