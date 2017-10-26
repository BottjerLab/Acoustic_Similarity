function [lasts, lastIdxs] = findLastborn(parentEv, childEv)
pIdxs = findParent(parentEv,childEv); %particular
lastIdxs = NaN(1,length(parentEv)); % last indices
for i = 1:length(parentEv)
    foo = find(pIdxs==i,1,'last');
    if ~isempty(foo)
        lastIdxs(i) = foo;
    end
end
lasts = childEv(lastIdxs(~isnan(lastIdxs))); %last approved syllables in a motif

end