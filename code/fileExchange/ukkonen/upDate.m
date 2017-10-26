function [s,k,leafcntr] = upDate(Txt, root, s, k, i, leafcntr)
%UPDATE Suffix Tree update function
%   Detailed explanation goes here

	% (s, (k, i-1)) is the canonical reference pair for the active point
	oldr = root;
	[endPoint,r] = test_and_split(Txt, s, k, i-1, Txt(i));
	
	while ~endPoint
		newState = State();		% create new leaf node
		newState.index=leafcntr;% increase leaf counter
		leafcntr=leafcntr+1;
		
		r=r.addTrans(Txt(i), i, length(Txt), newState);
		
		if oldr ~= root
			oldr.sLink = r;		% create new suffix link
		end
		
		oldr = r;
		
		[s,k] = canonize(Txt, s.sLink, k, i-1);
		[endPoint,r]= test_and_split(Txt, s, k, i-1, Txt(i));
	end
	
	if oldr ~= root
		oldr.sLink = s;
	end
end
