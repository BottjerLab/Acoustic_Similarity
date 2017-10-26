function stAr = initEmptyStructArray(fields, N)
% initialize struct array

if nargin < 2
    N = 0;
end
if isstr(fields),
    fields = {fields};
end
% fun with cell arrays
foo = cell(2,numel(fields));
[foo{1,:}] = fields{:};
[foo{2,:}] = deal(cell(N,1));
stAr = struct(foo{:});
end