function res = timeToString(tStamp) % in secconds
mins = floor(tStamp / 60);
secs = floor(mod(tStamp, 60));
leftover = mod(tStamp, 60) - secs;
res=sprintf('%02d:%05.2f', mins, secs + leftover);
end