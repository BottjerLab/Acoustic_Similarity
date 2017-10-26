function motifs = motif_search( T, str, depth, strlen, lower, upper, motifs )
%MOTIF_SEARCH search for potential motifs in a suffix tree
% T 	: root node of the suffix tree
% str	: original string from that the tree was built (incl. terminator)
% strlen: length of the substring / potential motif
% lower : lower bound for motif length
% upper : upper bound for motif lenght
% motifs: array to hold the motifs
	
	if ~isempty(T.transitions)
		len = length(T.transitions);

		for i=1:len,
			trans = T.transitions(i);
		
			if depth>1,
				if ~trans.state.isLeaf,
					% next node is not a leaf, recursive call
					strlen = strlen + (trans.right - trans.left + 1);
					motifs = motif_search(trans.state, str, depth-1, strlen, lower, upper, motifs);
					strlen = strlen - (trans.right - trans.left + 1);
				%else
				%	next node is a leaf -> motif occurs only once -> skip
				end
			end
			
			if ~trans.state.isLeaf && strlen+(trans.right-trans.left+1)>=lower,
				start=trans.left-strlen;
				stop =min(trans.right,trans.left+upper-strlen-1);
				motifs = [ motifs {str(start:stop) } ];
			end
		end
	end
	
	% cleanup the motifs that contain a terminator char
	mot=1;
	while mot<=length(motifs)
		if ~isempty(regexp(motifs{mot}, '[!§$%&/()=?]', 'once'))
			motifs(mot)=[];
		else
			mot=mot+1;
		end
	end
end