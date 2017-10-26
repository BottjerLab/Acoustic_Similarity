function res = versionNumber
% returns the version as a decimal number for compatibility purposes 
% (note: this doesn't work as well for 7.10 and up)'
    d=version;
    dots = strfind(d,'.');
    res = str2num(d(1:dots(2)-1));