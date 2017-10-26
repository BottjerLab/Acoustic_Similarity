function [s,k] = canonize(Txt, s, k, p)
%CANONIZE Suffix Tree canonize function
%   Detailed explanation goes here
	if p < k
		return		%return (s, k);
	end
	
	% find the t_k transition g'(s,(k',p'))=s' from s
	% k1 is k',  p1 is p'
	w1ands1 = s.getTrans(Txt(k));	% s --(w1)--> s1
	s1 = w1ands1.state;
	k1 = w1ands1.left;
	p1 = w1ands1.right;
	
	while p1-k1 <= p-k				% s --(w1)--> s1 ---> ...
		k = k + p1 - k1 + 1;		% remove |w1| chars from front of w
		s = s1;
		
		if k <= p
			w1ands1 = s.getTrans(Txt(k));		% s --(w1)--> s1
			s1 = w1ands1.state;
			k1 = w1ands1.left;
			p1 = w1ands1.right;
		end
	end
end
