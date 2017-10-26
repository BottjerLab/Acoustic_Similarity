function cm = interpColormap(colA, colB, N)
if nargin < 3,
    N = 128;
end
if nargin < 2 || isempty(colB)
    colB = [128 128 128]/255;
end
    
cm = zeros(N,3);
for ii = 1:3
    cm(:,ii) = linspace(colA(ii), colB(ii), N);
end
    