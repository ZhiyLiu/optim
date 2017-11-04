% Input:
%   - atomMatrix: cell array of m-rep primitives
%                or double array of m-rep primitives
%   - meanRadii: double array of mean radii 
%   - flag -> 1 : multiplication, flag -> 0 : division 
%   - weightMode : 1-meanRadius/2, 0-meanRadius
% Output:
%   - retAtomMatrix: cell array of m-rep primitives after adjustment on scales
% 
% !!!! size(atomMatrix) == size(meanRadii) !!!!
%


function retAtomMatrix = makeCommensurate(atomMatrix, meanRadii, flag, weightMode)

dims = (size(atomMatrix) == size(meanRadii));

if (~dims(1) || ~dims(2))
    error('    Error: the dimension of the first two parameters are not the same.');
end

[row, col] = size(atomMatrix);

for r = 1:row
    for c = 1:col
        retAtomMatrix{r, c} = adjustUnits(atomMatrix{r,c}, meanRadii(r,c), flag, weightMode);
    end
end


return;
