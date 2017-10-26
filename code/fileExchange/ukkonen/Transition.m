classdef Transition < handle
	%TRANSITION Suffix Tree Transition class
	%   A transition in the suffix tree is an edge that connects two states.
	%	It contains the starting character 'litera',
	%	label indices (left/right), and a reference to its child state.
	
	properties
		litera;
		left;
		right;
		isEmpty;
		state;
	end
	
	methods
		function obj=Transition(litera,left,right,state)
			obj.litera	= litera;
			obj.left	= left;
			obj.right	= right;
			obj.state	= state;
			obj.isEmpty	= right<left;
		end
	end
	
end

