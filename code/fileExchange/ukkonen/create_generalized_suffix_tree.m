function [root, sfstring] = create_generalized_suffix_tree(varargin)
% CREATE_GENERALIZED_SUFFIX_TREE produces a generalized suffix tree
%   for multiple strings using the CREATE_SUFFIX_TREE function
%	(Ukkonen '95)
%

	% string terminators to separate multiple strings
	terminators = '!§$%&/()=?';
	 
	% check how many strings have been submitted (currently 10 string at most)
	if nargin>length(terminators)
		disp(['can only handle ' num2str(length(terminators)) ' strings!']);
		disp('exiting!');
		return;
	end
	
	% concatenate the strings using these terminators
	sfstring = '';
	if nargin==1 && varargin{1}(end) == '!'
		sfstring = varargin{1};
	else
		for i=1:nargin,
			sfstring = [sfstring char(varargin{i}) char(terminators(i))];
		end
	end
	
	% create a suffix tree with the produced string
	root = create_suffix_tree(sfstring);
	
	% fix suffixes to produce the "real" generalized suffix tree output
	if nargin>1,
		fixSuffixes(root,sfstring);
	end
end

function fixSuffixes(T,Txt)
	len = length(T.transitions);
	if len>0		
		for i=1:len			
			trans=T.transitions(i);
			term = regexpi(Txt(trans.left:trans.right), '[!§$%&/()=?]');
			if ~isempty(term)
				trans.right=trans.left+term(1)-1;
			end
			fixSuffixes(T.transitions(i).state,Txt);
		end
	end
end
