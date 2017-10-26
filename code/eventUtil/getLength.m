function ret = getLength(region)
%GETLENGTH gets the duration of a region in seconds
%  
%   getClip(REGION) gets the length of a region (or array of regions) in seconds.  
ret = [region.stop] - [region.start];