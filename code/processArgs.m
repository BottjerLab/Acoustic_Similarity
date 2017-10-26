function currParams = processArgs(currParams, varargin)
%PROCESSARGS processes aliases and name/value pairs for parameters
%
%   PARAMS = processArgs(PARAMS, PARAMNAME, PARAMVALUE, ...) assigns the value
%   paramvalue to the corresponding paramname.  This is roughly equivalent to 
%   PARAMS.(PARAMNAME) = PARAMVALUE.  Any number of parameters can be
%   assigned this way, and new parameters can be defined as well.  
%   
%   If PARAMNAME refers to a member of a substructure, i.e.
%   'rough.Nfreqbands', then the parameter value will be assigned to that
%   substructure, i.e. params.rough.Nfreqbands.  In this way processArgs
%   gives additional functionality over .() assignment.
%   
%   PARAMS = processArgs(PARAMS, ALIASNAME, ...) expands to all the parameters
%   referred by the alias ALIASNAME, as defined in defineAliases.m.  

%   Aliases can be included in a list alongside parameter name/value pairs.
%   Later argument definitions take precedence.

if nargin == 1, return; end;
paramFields = fieldnames(currParams);

% TODO: if there are any aliases/flags, parse those first
varleft = varargin;

% define aliases
aliasDefs = defineAliases;

% first look for aliases
aliasnames = aliasDefs(:,1);
% check to make sure aliases don't collide with parameter names
badAliases = intersect(aliasnames, paramFields);

if ~isempty(badAliases),
    badAliases %#ok<NOPRT>
    error('processArgs:aliasNameReserved','Alias names collide with parameter names');
end

% do recursive expansion on aliases (aliases can contain other aliases)
ii = 1;
while ii <= numel(varleft)
    if ischar(varleft{ii})
        matchAlias = find(strcmp(varleft{ii},aliasnames));
        if ~isempty(matchAlias) 
            if numel(matchAlias) > 1, 
                error('processArgs:aliasMultiplyDefined', ...
                    'Alias ''%s'' multiply defined',aliasnames(matchAlias(1)));
            end
            expandedAlias = aliasDefs{matchAlias,2};
            varleft = [{varleft{1:(ii-1)}} expandedAlias {varleft{(ii+1):end}}];
        	ii = ii - 1;
        end
    end
    ii = ii + 1;
end
    
if mod(numel(varleft),2) == 1, error('Error: Unbound alias or unpaired argument requires another value.'); end;

% now process normal pairs
argnames = varleft(1:2:end);
argvals = varleft(2:2:end);
for ii = 1:numel(argnames)
    iArg = argnames{ii};
    iVal = argvals{ii};
    % is this assignable to a substruct?
    if isempty(strfind(iArg,'.')) 
        if isempty(strfind(paramFields,iArg))
            error('Argument %s does not appear in list',iArg);
        else
            currParams.(iArg) = iVal;
        end
    else
        % recursively unravel this argument...
        [rootStruct, leafField] = strtok(iArg,'.');
        leafField(1) = '';
        currParams.(rootStruct) = ...
            processArgs(currParams.(rootStruct), leafField, iVal);
    end
end

end