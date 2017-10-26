%   The algorithms implemented by Alexander Vezhnevets aka Vezhnick
%   <a>href="mailto:vezhnick@gmail.com">vezhnick@gmail.com</a>
%
%   Copyright (C) 2005, Vezhnevets Alexander
%   vezhnick@gmail.com
%   
%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.
%
%   tree_node_w Implements the constructor for tree_node_w class, that
%   imlements classification tree
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    tree_node = tree_node_w(max_split)
%    ---------------------------------------------------------------------------------
%    Arguments:
%           max_split - maximum number of splits in the tree
%    Return:
%           tree_node - object of tree_node_w class

% 2012 Stathis Fotiadis, no (C)
% Changes to object saving and loading functions

classdef tree_node_w %< handle
    properties %(SetAccess = public)
        left_constrain
        right_constrain
        dim
        max_split
        parent
    end
    
    methods 
        function tn = tree_node_w(ms)
            tn.max_split = ms;
            %tree_node = class(tree_node, 'tree_node_w') ;
        end

         function A = saveobj(obj)
           A.left_constrain  = obj.left_constrain;
           A.right_constrain = obj.right_constrain;
           A.dim             = obj.dim;
           A.max_split       = obj.max_split;
           A.parent          = obj.parent;
         end 
    end
    
    methods (Static)
      function obj = loadobj(A)
         obj = tree_node_w(A.max_split);
         obj.left_constrain  = A.left_constrain;
         obj.right_constrain = A.right_constrain;
         obj.dim             = A.dim;
         obj.parent          = A.parent;
      end
    end
    
end