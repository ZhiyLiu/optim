% Input:
%   - atomMatrix: cell array of m-rep primitives
%
% Output:
%   - radii: double array of radius
% 
%


function radii = getRadii(atomMatrix)

[row, col] = size(atomMatrix);

for r = 1:row
    for c = 1:col
        radii(r, c) = get(atomMatrix{r,c}, 'r');
    end
end


return;