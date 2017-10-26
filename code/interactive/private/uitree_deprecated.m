function [tree, container] = uitree_deprecated(varargin)
% This function is undocumented and will change in a future release

%   Copyright 2003-2012 The MathWorks, Inc.

%   Release: R14. This feature will not work in previous versions of MATLAB.

% Warn that this old code path is no longer supported.
warn = warning('query', 'MATLAB:uitree:DeprecatedFunction');
if isequal(warn.state, 'on')
    warning(message('MATLAB:uitree:DeprecatedFunction'));
end

%% Setup and P-V parsing.

error(javachk('jvm'));
error(nargoutchk(0, 2, nargout));

fig = [];
numargs = nargin;

if (nargin > 0 && isscalar(varargin{1}) && ishandle(varargin{1}))
    if ~ishghandle(handle(varargin{1}), 'figure')
        error(message('MATLAB:uitree:InvalidFigureHandle'));
    end
    fig = varargin{1};
    varargin = varargin(2:end);
    numargs = numargs - 1;
end

% RootFound = false;
root   = [];
expfcn = [];
selfcn = [];
pos    = [];
% parent = [];

if (numargs == 1)
    error(message('MATLAB:uitree:InvalidNumInputs'));
end

for i = 1:2:numargs-1
    if ~ischar(varargin{i})
        error(message('MATLAB:uitree:UnrecognizedParameter'));

    end
    switch lower(varargin{i})
        case 'root'
            root = varargin{i+1};
        case 'expandfcn'
            expfcn = varargin{i+1};
        case 'selectionchangefcn'
            selfcn = varargin{i+1};
        case 'parent'
            if ishandle(varargin{i+1})
                f = varargin{i+1};
                if ishghandle(handle(f), 'figure')
                    fig = f;
                end
            end
        case 'position'
            p = varargin{i+1};
            if isnumeric(p) && (length(p) == 4)
                pos = p;
            end
        otherwise
            error('MATLAB:uitree:UnknownParameter', ['Unrecognized parameter: ', varargin{i}]);
    end
end

if isempty(expfcn)
    [root, expfcn] = processNode(root);
else
    root = processNode(root);
end

% This should not be tagged for EDT invocation. This is a FigureChild.
% Its public methods expect main thread invocation which subsequently post
% stuff to EDT.
tree_h = com.mathworks.hg.peer.UITreePeer; 
tree_h.setRoot(root);

if isempty(fig)
    fig = gcf;
end

if isempty(pos)
    figpos = get(fig, 'Position');
    pos =  [0 0 min(figpos(3), 200) figpos(4)];
end

[tree, container] = javacomponentfigurechild_helper(tree_h,pos, fig);

if ~isempty(expfcn)
    set(tree, 'NodeExpandedCallback', {@nodeExpanded, tree, expfcn});
end

if ~isempty(selfcn)
    set(tree, 'NodeSelectedCallback', {@nodeSelected, tree, selfcn});
end

end

%% -----------------------------------------------------
function nodeExpanded(src, evd, tree, expfcn)                           %#ok

% tree = handle(src);
% evdsrc = evd.getSource;

evdnode  = evd.getCurrentNode;
% indices = [];

if ~tree.isLoaded(evdnode)
    value = evdnode.getValue;

    % <call a user function(value) which returns uitreenodes>;
    cbk = expfcn;
    if iscell(cbk)
        childnodes = feval(cbk{1}, tree, value, cbk{2:end});
    else
        childnodes = feval(cbk, tree, value);
    end

    if (length(childnodes) == 1)
        % Then we dont have an array of nodes. Create an array.
        chnodes = childnodes;
        childnodes = javaArray('com.mathworks.hg.peer.UITreeNode', 1);
        childnodes(1) = java(chnodes);
    end

    tree.add(evdnode, childnodes);
    tree.setLoaded(evdnode, true);
end

end

%% -----------------------------------------------------
function nodeSelected(src, evd, tree, selfcn)                           %#ok
cbk = selfcn;
hgfeval(cbk, tree, evd);

end

%% -----------------------------------------------------
function [node, expfcn] = processNode(root)
expfcn = [];

if isempty(root) || isa(root, 'com.mathworks.hg.peer.UITreeNode') || ...
        isa(root, 'javahandle.com.mathworks.hg.peer.UITreeNode')
    node = root;
elseif ishghandle(root)
    % Try to process as an HG object.
    try
        %In the HG tree rendering, the java side is intelligent enough to
        %pick an icon based on the class of the object that we are sending
        %in. It is absolutely essential that we send in the handle to java
        %because the icon picking is based on the bean adapter's display
        %name. 
        node = uitreenode_deprecated(handle(root), get(root, 'Type'), ...
            [], isempty(get(0, 'Children')));
    catch ex %#ok mlint
        node = [];
    end
    expfcn = @hgBrowser;
elseif ismodel(root)
    % Try to process as an open Simulink system

    % TODO if there is an open simulink system and a directory on the path with
    % the same name, the system will hide the directory. Perhaps we should
    % warn about this.
    try
        h = handle(get_param(root,'Handle'));
        % TODO we pass root to the tree as a string,
        % it would be better if we could just pass the
        % handle up
        node = uitreenode_deprecated(root, get(h, 'Name'), ...
            [], isempty(h.getHierarchicalChildren));
    catch ex %#ok mlint
        node = [];
    end
    expfcn = @mdlBrowser;
elseif ischar(root)
    % Try to process this as a directory structure.
    try
        iconpath = [matlabroot, '/toolbox/matlab/icons/foldericon.gif'];
        node = uitreenode_deprecated(root, root, iconpath, ~isdir(root));
    catch ex %#ok mlint
        node = [];
    end
    expfcn = @dirBrowser;
else
    node = [];
end

end

%% -----------------------------------------------------
function nodes = hgBrowser(tree, value)                                 %#ok

try
    count = 0;
    parent = handle(value);
    ch = parent.children;

    for i=1:length(ch)
        count = count+1;
        nodes(count) = uitreenode_deprecated(handle(ch(i)), get(ch(i), 'Type'), [], ...
            isempty(get(ch(i), 'Children')));
    end
catch
    error(message('MATLAB:uitree:UnknownNodeType'));
end

if (count == 0)
    nodes = [];
end

end

%% -----------------------------------------------------
function nodes = mdlBrowser(tree, value)                                %#ok

try
    count = 0;
    parent = handle(get_param(value,'Handle'));
    ch = parent.getHierarchicalChildren;

    for i=1:length(ch)
        if isempty(findstr(class(ch(i)),'SubSystem'))
            % not a subsystem
        else
            % is a subsystem
            count = count+1;
            descr = get(ch(i),'Name');
            isleaf = true;
            cch =  ch(i).getHierarchicalChildren;
            if ~isempty(cch)
                for j = 1:length(cch)
                    if ~isempty(findstr(class(cch(j)),'SubSystem'))
                        isleaf = false;
                        break;
                    end
                end
            end
            nodes(count) = uitreenode_deprecated([value '/' descr], descr, [], ...
                isleaf);
        end
    end
catch
    error(message('MATLAB:uitree:UnknownNodeType'));
end

if (count == 0)
    nodes = [];
end

end


%% -----------------------------------------------------
function nodes = dirBrowser(tree, value)                                %#ok

try
    count = 0;
    ch = dir(value);

    for i=1:length(ch)
        if (any(strcmp(ch(i).name, {'.', '..', ''})) == 0)
            count = count + 1;
            if ch(i).isdir
                iconpath = [matlabroot, '/toolbox/matlab/icons/foldericon.gif'];
            else
                iconpath = [matlabroot, '/toolbox/matlab/icons/pageicon.gif'];
            end
            nodes(count) = uitreenode_deprecated([value, ch(i).name, filesep], ...
                ch(i).name, iconpath, ~ch(i).isdir);
        end
    end
catch
    error(message('MATLAB:uitree:UnknownNodeType'));
end

if (count == 0)
    nodes = [];
end

end

%% -----------------------------------------------------
function yesno = ismodel(input)
yesno = false;

try
    if is_simulink_loaded 
        get_param(input,'handle');
        yesno = true;
    end
end

end
