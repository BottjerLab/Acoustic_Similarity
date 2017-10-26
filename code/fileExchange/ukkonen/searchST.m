function indices = searchST(T,Str,Qry)
%searchST checks if a string is in a suffix tree and returns occurrence indices
%   suffix tree root: T
%	original string : Str
%   query string    : Qry
	
	if ~isempty(T.transitions)
		curr_state = T;			% store current state for the iterative approach
		qrylen = length(Qry);	% store query length
		indices=[];				% initialize the indices array
		i = 1;					% position variable
		
		while i<=qrylen
			% get array index of next transition by checking the character
			if ~curr_state.isLeaf,
				pos = regexp(char(curr_state.transitions.litera)', Qry(i), 'once');
			else
				pos = [];
			end
			
			% if( state.transition(X).litera != Qry(i) ) -> break -> no match
			if isempty(pos),
				break;
			end
			
			trans  = curr_state.transitions(pos);	% get current transition
			strlen = trans.right-trans.left+1;		% and its label length
			
			if length(Qry(i:end))<=strlen
				% compare the substring and the query
				result = strcmp( Qry(i:end), Str( trans.left:trans.left+length(Qry(i:end))-1) );
				
				% if query was found -> get the occurrence indices
				if result
					if trans.state.isLeaf
						% next state is a leaf -> only one occurrence
						% compute the occurrence index
						indices = trans.left + length(Qry(i:end)) - qrylen;
					else
						% next state is not a leaf -> multiple occurrences
						% need to traverse the subtree:
						indices = getIndices( trans.state );
					end
				end
				
				break;	% break the loop
			else
				% query string longer than edge string:
				% update the current state for next iteration step
				curr_state = trans.state;
				% update the index i
				i = i+strlen;
			end
		end
	end
end

% traverse the subtree (from state "state") in order to compute occurrence indices
function indices = getIndices(state)
	indices = [];
		
	for j=1:length(state.transitions)
		if state.transitions(j).state.isLeaf
			indices = [ indices state.transitions(j).state.index ];
		else
			indices = [ indices getIndices(state.transitions(j).state) ];
		end
	end
end


