% function sumM = mrepSum(M1, deltaM2)
% 
% Compute differences of corresponding atoms in deltaM2 and M1
% !!!!size(deltaM2) == size(M1)!!!
%
% Input:
%   - deltaM2, M1: cell array of m-rep primitives (either tube or quad)
% Output:
%   - sumM: cell array (sumM = M1 + deltaM2)
%           NOTE: NOT commutative!!!!
% 

function sumM = mrepSum(M1, deltaM2)

[row, col] = size(M1);
[row1, col1] = size(deltaM2);
if (col ~= col1 | row ~= row1) 
    error('    Error: The dimension of two input cell arrays are not the same.');
end

for c =1:col
    for r =1:row
        sumM{r,c} = M1{r,c} + deltaM2{r,c};
    end
end

return;
