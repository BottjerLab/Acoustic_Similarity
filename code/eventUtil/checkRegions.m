function [verified, repaired] = checkRegions(regions,fs)

% tolerate off by one errors - nobody's perfect
verified = abs([regions.start] * fs -      [regions.idxStart]) <= 1 & ...
           abs([regions.stop ] * fs - ceil([regions.idxStop] )) <=1;

if any(~verified) % use the timestamps as veridical
    repaired = eventFromTimes([regions.start],[regions.stop],fs);
    
    % bring the old labels to the new structure
    foo = {regions.type}; 
    [repaired.type] = foo{:};
else
    repaired = regions;
end
