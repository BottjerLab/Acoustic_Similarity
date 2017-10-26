function [subStrSorted, countsSorted, locations] = mostCommonSubstring(string,N, M)
% returns the most common substrings of length N, with more than M
% occurrences
    if nargin < 3
        M = 1;
    end
    [subStr, counts, locations] = n_gram(string, N);
    [countsSorted, sortIdx] = sort(counts, 'descend');
    subStrSorted = subStr(sortIdx);
    subStrSorted = subStrSorted(countsSorted > M)';
    countsSorted = countsSorted(countsSorted > M)';
    rIdx = zeros(1,numel(sortIdx)); 
    rIdx(sortIdx) = 1:numel(sortIdx);
    locations = rIdx(locations);
    
end

function [subStrings, counts, index] = n_gram(fullString, N)
  if (N == 1)
    [subStrings, rIdx, index] = unique(cellstr(fullString.'));  %.'# Simple case
    subStrings{cellfun('isempty',subStrings)} = ' ';
  else
    nString = numel(fullString);
    index = hankel(1:(nString-N+1), (nString-N+1):nString);
    [subStrings, rIdx, index] = unique(cellstr(fullString(index)));
    
    % make sure substrings have trailing spaces
    for ii = 1:numel(subStrings)
        if numel(subStrings{ii}) ~= N
            subStrings{ii} = [subStrings{ii} ' ']; %assume only single spaces exist
        end
    end
  end
  counts = accumarray(index, 1);
end
