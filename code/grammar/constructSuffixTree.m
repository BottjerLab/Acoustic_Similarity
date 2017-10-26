function treeArray = constructSuffixTree(sequence)
% initialize variables
active = struct('node',1,'edge',[],'edgechar',[],'length',0); % node = index in treeArray
% edge is the index inside the node's children of the edge
remainder = 0;

% first node is root node
% here the # 0 in bounds represents the end
treeArray = struct('bounds',{[0,0], [1,0]},'children',{2, []},'link',[]);
freshNode = struct('bounds',[],'children',[],'link',[]);
rootNodeIdx = 1;

for currPos = 2:numel(sequence)
    nodeCreated = false;
    
    % is there another node that begins with the same character?
    addNodes = false;
    if active.length == 0 % we're at the root or end of a node
        % is there a new edge which contains our character
        addNodes = isempty(findEdgeInit(treeArray(active.node), sequence(currPos)));
    else
        % is the next character in the edge this character?
        nodeInTree = treeArray(treeArray(active.node).children(active.edge));
        nextCharOfEdge = sequence(nodeInTree.bounds(1) + active.length + 1);
        addNodes = (sequence(currPos) ~= nextCharOfEdge);
    end
    
    if ~addNodes % this suffix is an extension of an existing prefix
        active.length = active.length + 1;
        remainder = remainder + 1;
        nodeInTree = treeArray(treeArray(active.node).children(active.edge));
        if active.length == lenOf(nodeInTree);
            % if we've run off the current edge in terms of characters, we have
            % to move to that node
            active.node = active.edge;
            active.edge = [];
            active.edgechar = [];
            active.length = 0;
        end        
    else % we need to add at least one node
        while remainder > 0
            % this is a new character, not in the start of this node
            if active.node == rootNodeIdx && ~isempty(active.edge) 

                % <--------------------------
                % add new child edge/node
                treeArray(end+1) = freshNode;
                treeArray(end).bounds = [currPos -1];
                % add this node to the list of the root's children
                treeArray(rootNodeIdx).children(end+1) = numel(treeArray);
                remainder = remainder - 1;
            else % this is a new element that we need to split an edge off of
                % split the edges - prepare new nodes
                nNodes = length(treeArray);
                newNode1 = freshNode;
                newNode1.bounds = [currPos -1];
                newNode2 = freshNode;
                % insert new nodes
                newNode2.bounds = [treeArray(active.node).bounds(1) + active.length, -1];
                treeArray(end+1:end+2) = [newNode1, newNode2];
                remainder = remainder - 1;
                
                % trim the edge that is split
                splitNodeIdx = treeArray(active.node).children(active.edge);
                treeArray(splitNodeIdx).bounds(2) = ...
                    treeArray(active.node).bounds(1) + active.length - 1;
                % connect edges
                treeArray(splitNodeIdx).children = [treeArray(active.node).children nNodes + 1 nNodes + 2];
                
                % update active locations
                
                if active.node == rootNodeIdx
                    newChar = (currPos - remainder + 1);
                    active.length = active.length - 1;
                    active.edge = findEdgeInit(treeArray(rootNodeIdx));
                    
                else
                end
            end
        end
    end
    
end

    function ret = lenOf(myNode)
        ret = diff(myNode.length);
        if ret < 0; ret = ret + currPos + 1;end
    end

    function ret = findEdgeInit(treeNode, char)
        % returns [] if not found
        ret = [];
        for ll = 1:numel(treeNode.children)
            if sequence(treeArray(treeNode.children(ii)).bounds(1)) == char
                ret = ll; return;
            end
        end
    end

end