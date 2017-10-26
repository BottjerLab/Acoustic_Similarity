function root = create_suffix_tree( Txt )
%CREATE_SUFFIX_TREE create a suffix tree (Ukkonen '95)
%   from E.Ukkonen, On-Line Construction of Suffix Trees ***
%                 Algorithmica 14(3) pp 249-260, 1995  ***
%					(U. Helsinki, Finland)
	
% 	disp(['Creating suffix tree (string length: ' ...
% 		num2str(length(Txt)) '). This might take a while!']);
	
	root = State();
	bt = State();			% bt (bottom or _|_)
	
	% create transitions for all possible chars from bt to root
	for i=1:length(Txt)
		bt = bt.addTrans(Txt(i),i,i,root);
	end
	
	root.sLink = bt;	% create suffix link
	s=root;
	k=1;
	
	leafcntr = 1;
	for i=1:length(Txt)
		if mod(i,1000)==0,
			disp(['  `- processing ' num2str(i) '-th symbol']);
		end
		
		% Txt string (and the root node handler) are needed for the methods
		[s,k,leafcntr] = upDate(Txt, root, s, k, i, leafcntr);	% (s,k) < - upDate(...)
		[s,k] = canonize(Txt, s, k, i);							% (s,k) < - canonize(...)
	end
	
% 	disp('Suffix tree has been created! Thanks for waiting.');
end
