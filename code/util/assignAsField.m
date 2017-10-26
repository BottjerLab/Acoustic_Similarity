function structArray = assignAsField(structArray, fieldName, fieldData)
assert(numel(fieldData) == numel(structArray));
num2cell(fieldData);
[structArray.(fieldName)] = ans{:};