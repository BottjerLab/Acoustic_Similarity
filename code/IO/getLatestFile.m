function [fil, filExist] = getLatestFile(glob)
fileDir = fileparts(glob);
results = dir(glob); 
fil = []; filExist = false;
if ~isempty(results)
    results = sortBy(results,'datenum', 'descend');
    fil = [fileDir filesep results(1).name];
    filExist = true;
end
end