% 
% srepQuadPrims is ONE or a ROW or COLUMN vector of QuadSymPrimitives
%
% Output srepQuadVec :  (num_atoms_srepQuadPrims) x (dim of atom 9 or 11)

% !! scan Row-wise

function srepQuadVec = cVec(srepQuadPrims)

[r, c] = size(srepQuadPrims);
if ~(r == 1 || c == 1)
    error('Input parameter must be either row or column vector');
    % Why?
end

srepQuadVec = [];

for i = 1:r
    for j = 1:c
        srepQuadPrim = srepQuadPrims(i,j);
        v = [srepQuadPrim.pos' srepQuadPrim.r  reshape(srepQuadPrim.U, 1, [])];
        srepQuadVec = [srepQuadVec; v];
    end
end
