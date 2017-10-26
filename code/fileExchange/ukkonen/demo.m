%% Demo file for SUFFIX TREE implementation
%
% Ukkonen 1995: http://www.springerlink.com/content/kq55005qu6479276/
% source-code inspiration: http://www.allisons.org/ll/AlgDS/Tree/Suffix/
%

% your original string. make sure it is terminated with a uniqe char!
%txt = 'mississippi$'
%txt = 'hattivattidattibatti$'
%txt = 'bananas$'
%txt = 'ababbb$'
txt = 'mississippi$';

% create the suffix tree
[root] = create_suffix_tree(txt);

% display the results
showST(root, txt, 0);



%% create generalized suffix tree and display it
[gst_root,treestring] = create_generalized_suffix_tree('abab','baba','abba');

showST(gst_root, treestring, 1);



%% check if a string pattern occurrs in the original string / the tree 
query = 'issi';
occ_ind = searchST(root,txt,query);

disp( ['Query string: "' query '"'] );
disp( ['No. of occurences: ' num2str(length(occ_ind)) ] );
disp( ['Corresponding indices: ' num2str(occ_ind)]);


