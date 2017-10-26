function [endPoint,state] = test_and_split(Txt, s, k, p, t)
%TEST_AND_SPLIT Suffix Tree test-and-split function
%   Detailed explanation goes here
	
	if k<=p
		% find the t_k transition g'(s,(k',p'))=s' from s
		% k1 is k'  p1 is p'
		w1ands1 = s.getTrans(Txt(k));	% s --(w1)--> s1
		s1 = w1ands1.state;				% extract	the state
		k1 = w1ands1.left;				% ...		the left index
		p1 = w1ands1.right;				% ...		the right index
			
		if t==Txt(k1 + p - k + 1)
			% return the pair(true, s);
			endPoint = true;
			state = s;
		else
			% replace the t_k transition (above) by following two transitions:
			% s --trans1--> r --trans2--> s1
			r = State();

 			% replace the existing transition first	(s ---> r)
			%w1ands1.left = k1; w1ands1.right = k1+p-k; w1ands1.state = r;
			s = s.addTrans(Txt(k1), k1, k1+p-k, r);
			
			% and add another one (r ---> s1)
			r = r.addTrans(Txt(k1+p-k+1), k1+p-k+1, p1, s1);
			
			% return the pair(false, r)
			endPoint = false;
			state = r;
		end
	else
		% is there a t-transition from s ?
		%return the pair(s[t] != null, s);
		endPoint = s.getTrans(t) ~= 0;
		state = s;
	end
end
