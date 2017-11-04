% ==============================================================
% Function valStr = readStringLines(keyStr, modelArray)
%
%   Input
%   modelArray: a character array of m3d model file
%   keyStr: a key string with which lines that we want to
%               find in modelArray start
%   patternStr: A pattern to match on the string returned.
%   Output: values in lines corresponding to the key as a cell array of
%   strings
% --------------------------------------------------------------
function valStr = readStringLines(keyStr, modelArray, patternStr)
    Lines = modelArray(strmatch(keyStr, modelArray), :);
    [nRow, nCol] = size(Lines);
    valStr = {};
    for i = 1:nRow
        newVal = sscanf(Lines(i,:), [keyStr, ' = %s;\n']);
        pos = strfind(newVal, patternStr);
        if(size(pos,1) ~= 0 )
            % pattern match successful.
            % Remove the ending semicolon :
            newVal  = newVal(1:end-1);
            valStr = { valStr{:} newVal };
        end
    end
return;
