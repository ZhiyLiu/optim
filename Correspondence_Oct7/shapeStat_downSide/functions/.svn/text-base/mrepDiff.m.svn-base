% function diffM = mrepDiff(M2, M1)
% 
% Compute differences of corresponding atoms in M2 and M1
% !!!!size(M2) == size(M1)!!!
%
% Input:
%   - M2, M1: cell array of m-rep primitives (either tube or quad)
% Output:
%   - diffM: cell array (diffM = M2 - M1)
% 

function diffM = mrepDiff(M2, M1)

[row, col] = size(M1);
[row1, col1] = size(M2);
if (col ~= col1 | row ~= row1) 
    disp('    Error: The dimension of two input cell arrays are not the same.');
    disp('');
    return;
end

for c =1:col
    for r =1:row
        diffM{r,c} = M2{r,c} - M1{r,c};
    end
end

return;
