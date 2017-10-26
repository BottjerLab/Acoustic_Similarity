function clearSummaries
if strcmp(lower(input('Are you sure you want to clear summaries? ', 's')), 'y')
    delete('summaries\contents*.mat');
end