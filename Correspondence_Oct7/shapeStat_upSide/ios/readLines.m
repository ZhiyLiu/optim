% ==============================================================
% Function valStr = readLines(keyStr, modelArray)
%
%   Input
%   modelArray: a character array of m3d model file
%   Output: values in lines corresponding to the key
% --------------------------------------------------------------
function valStr = readLines(keyStr, modelArray)
    Lines = modelArray(strmatch(keyStr, modelArray), :);
    [nRow, nCol] = size(Lines);
    valStr = [];
    for i = 1:nRow
        newVal = sscanf(Lines(i,:), [keyStr, ' = %f;\n']);
        valStr = [valStr; newVal];
    end
return;
