function warnfiguredialog(functionName)

%   Copyright 2006-2012 The MathWorks, Inc.

% Shows the warning in -nodisplay and -noFigureWindows mode
% Shows the error in -nojvm mode

% (usejava('jvm') will return true in -nodisplay and -noFigureWindows mode,
% but false in -nojvm mode)

if ~usejava('jvm')
    error(message('MATLAB:HandleGraphics:noJVM'));
end

% @TODO When we are ready to deprecate dialogs in NoFigureWindows modes
% other than nojvm, use this check to throw errors:
% if feature('NoFigureWindows')

% if ~usejava('awt') % adding warnings for the noawt
if ~feature('ShowFigureWindows')
    warningId = strcat('MATLAB:', functionName, ':DeprecatedBehavior');

	% @HACK HACK alert: if you are modifying this message, please also modify
	% the same message in the following files:
	%  * matlab/toolbox/matlab/graphics/private/warnfiguredialog.m
    warning(warningId, sprintf(getString(message('MATLAB:warnfiguredialog:DeprecatedBehavior'))));
end
