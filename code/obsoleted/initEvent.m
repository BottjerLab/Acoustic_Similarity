function events = initEvent
% flag: isn't that helpful, should be replaced/obsoleted
% initialize sparse coding of status, using events, which will be a struct
% array
fields = {'type','start','stop','idxStart','idxStop'};

% fun with cell arrays
foo = cell(2,numel(fields));
[foo{1,:}] = fields{:};
[foo{2,:}] = deal(NaN);
events = struct(foo{:});
end