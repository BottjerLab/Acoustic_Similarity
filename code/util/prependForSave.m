function ret = prependForSave(prefix, orig)

% prepends prefix only to the filename portion of the file, not the
% directory
[dired, fil, ext] = fileparts(orig);
fil = [prefix fil];
ret = fullfile(dired,[fil,ext]);
end