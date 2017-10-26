function dst = copyPartialStruct(src, dst, fieldsToCopy)
% warning: only works for numeric fields, single structures at at time

for kk = 1:numel(fieldsToCopy)
    foo = fieldsToCopy{kk};
    num2cell([src.(foo)]);
    [dst.(foo)]= ans{:};
end