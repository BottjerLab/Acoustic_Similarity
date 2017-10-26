function [structArray, structOfArrays] = readDataFile(dataFile)
% READDATATFILE read space/tab delimited files

% see RECORDINGDATA.txt for an example of acceptable files
% todo: see if we can handle xls files, by some kind of underhanded
% conversion

if nargin < 1
    [dataFile dataPath] = uigetfile('Please select data file to read');
    dataFile = [dataPath dataFile];
end

if ~exist(dataFile,'file')
    error('File does not exist!');
end
fid = fopen(dataFile);

% read the header if there is one
header = fgetl(fid); % TODO: check if empty?
fields = [];
if header(1) == '%'
    fields = regexp(header,'<(\w+)>','tokens');
    fields = [fields{:}];
    nFields = numel(fields);
else
    nFields = 0; 
    while ~isspace(header)
        [~,header]=strtok(header,' '); 
        nFields = nFields + 1;
    end
    fields = strcat('Field', num2str([1:nFields]'));
    warning('No header, reading as strings');
    fseek(fid,0,-1);
end

forstr = repmat('%s ', 1, nFields); 
forstr(end) = '';

cellArray = textscan(fid,forstr,'CommentStyle','%');
nR = numel(cellArray{1});

structArray = initEmptyStructArray(fields, nR);
structOfArrays = initEmptyStructArray(fields, 1);
for ii = 1:nFields
    dataCell = cellArray{ii};
    % is it numeric?
    numericTest = str2double(dataCell);
    if all(~isnan(numericTest)) % yes, it is
        dataCell = num2cell(numericTest); 
        structOfArrays.(fields{ii}) = numericTest;
    else
        structOfArrays.(fields{ii}) = dataCell;
    end
    [structArray.(fields{ii})] = dataCell{:};
end
fclose(fid);
end