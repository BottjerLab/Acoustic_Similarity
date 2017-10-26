function depth=showST(T,Txt,indent)
%SHOWST Show a suffix tree
%   takes the suffix tree T
%	the original string   Txt
%	and the counter	      indent (original value to be 0)
	
	len = length(T.transitions);	% number of transitions T has
	
	if len>0
		spaces = '';
		if indent==0
			disp(char(zeros(1,length(Txt)+2)+'-'));
			disp([' ' Txt ' ']);
			disp(char(zeros(1,length(Txt)+2)+'-'));
		else
			spaces = char(zeros(1,indent*4));
		end
		
		for i=1:len			
			trans=T.transitions(i);
			if ~trans.state.isLeaf,
				disp( [spaces '  |--(' num2str(trans.left) ':' Txt(trans.left:trans.right) ') '  ] );
				indent = indent + 1;
				showST(T.transitions(i).state,Txt,indent);
				indent = indent - 1;
			else
				disp( [spaces '  |--(' num2str(trans.left) ':' Txt(trans.left:trans.right) ') [' num2str(trans.state.index) ']'  ] );
			end
		end		
	end
end
