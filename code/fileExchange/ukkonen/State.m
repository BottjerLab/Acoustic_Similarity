classdef State < handle
	%STATE Suffix Tree State class
	%   A suffix tree state contains the isLeaf flag, a leaf node index,
	%	a set of outgoing transitions (edges), and a suffix link.
	properties
		isLeaf;		sLink; 		transitions;	index;
	end
	
	methods
		function obj=State()
			obj.isLeaf=1;
			obj.transitions=[];
		end
		
		% add a transition to the state % this['a'] >--(left..right)--> s
		function obj=addTrans(obj,litera,left,right,state)
			% is there a transitions already stored?
			if obj.isLeaf == 0
				% check if there is a transition with the corresponding 'litera'
				[foo,loc]=ismember(litera,{obj.transitions.litera});			
				if loc>0
					% 'litera' transition available, overwrite the values
					obj.transitions(loc).left=left;
					obj.transitions(loc).right=right;
					obj.transitions(loc).state=state;
				else
					% store the new 'litera' transition 
					obj.transitions=[obj.transitions; Transition(litera,left,right,state)];
				end
			else
				% no transitions yet stored, store this one now
				obj.transitions=[obj.transitions; Transition(litera,left,right,state)];
				obj.isLeaf=0;
			end
		end
		
		% get the transition to a given litera... trans = s[Txt(k)]
		function trans=getTrans(obj,litera)
			if obj.isLeaf	% if there are no transitions at all
				trans=0;
				return;
			end
			
			% check for the transition with the char 'Txt(k)'
			[foo,loc]=ismember(litera,{obj.transitions.litera});			
			if loc>0
				trans=obj.transitions(loc);		% transition is available
			else
				trans=0;						% no transition stored
			end				
		end
	end
end
