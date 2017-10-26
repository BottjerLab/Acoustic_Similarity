function nd = minEditDistance(a, b)
if length(a) < length(b)
    nd = minEditDistance(b,a);
    return;
end

% a is longer than b
nOffsets = length(a) - length(b) + 1;
allNd = zeros(1,nOffsets);
for ii = 1:nOffsets
    allNd(ii) = sum(xor(a(ii:ii+length(b)-1),b));
end
nd = min(allNd)/length(b);