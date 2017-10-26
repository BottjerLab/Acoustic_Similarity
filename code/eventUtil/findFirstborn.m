function [firsts, firstIdxs] = findFirstborn(parentEv, childEv)
pIdxs = findParent(parentEv,childEv); %particular
firstIdxs = NaN(1,length(parentEv)); % first indices
for i = 1:length(parentEv)
    foo = find(pIdxs==i,1,'first');
    if ~isempty(foo)
        firstIdxs(i) = foo;
    end
end
firsts = childEv(firstIdxs(~isnan(firstIdxs))); %first approved syllables in a motif

end