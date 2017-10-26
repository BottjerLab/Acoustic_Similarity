function h = uitab_deprecated(varargin)
% This function is undocumented and will change in a future release

%   Copyright 2004-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/03/09 07:07:42 $

%   Release: R14SP2. This feature will not work in MATLAB R13 and before.

% Warn that this old code path is no longer supported.
warn = warning('query', 'MATLAB:uitab:DeprecatedFunction');
if isequal(warn.state, 'on')
    warning(message('MATLAB:uitab:DeprecatedFunction'));
end

error(javachk('swing'));

h = uitools.uitab(varargin{:});
h = double(h);
