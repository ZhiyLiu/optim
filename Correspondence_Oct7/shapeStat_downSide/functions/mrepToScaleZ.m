% Input:
%   - atomMatrix: cell array of m-rep primitives
%                 or double array of m-rep primitives
%   - weScaleVec: column of float scales for z dimension per each model
%   
% Output:
%   - retAtomMatrix: cell array of m-rep primitives after adjustment 
% 
%   !!!! size(atomMatrix, 1) = size(weScaleVec, 1) !!!!!
%


function retAtomMatrix = mrepToScaleZ(atomMatrix, weScaleVec)


[row, col] = size(atomMatrix);

if ( row ~= length(weScaleVec))
    disp('    Error: the number of rows of first and the second parameters are not the same.');
    disp(' ');
    return;
end



for r = 1:row
    for c = 1:col
        retAtomMatrix{r, c} = toScaleZ(atomMatrix{r,c}, weScaleVec(r));
    end
end


return;