function boolVec = stateVector(ev, subEv, resol)
if nargin < 3
    resol = 0.001;
end
ii = 1;
nTimeSteps = floor((ev.stop - ev.start)/resol);
boolVec = false(1, floor((ev.stop - ev.start)/resol));
for tt = ev.start:resol:ev.stop
    boolVec(ii) = any([subEv.start] <= tt & [subEv.stop] >= tt);
    ii = ii+1;
end
end