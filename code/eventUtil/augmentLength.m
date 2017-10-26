function regions = augmentLength(regions)
if isempty(regions) || isfield(regions,'length'), return; end;

foo = num2cell([regions.stop] - [regions.start]);
[regions.length] =  foo{:};