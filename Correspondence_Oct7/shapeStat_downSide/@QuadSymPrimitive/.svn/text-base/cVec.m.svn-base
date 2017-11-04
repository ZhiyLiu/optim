% 
% quadSymPrims is ONE or a ROW or COLUMN vector of QuadSymPrimitives
%
% Output quadSymVec :  (num_atoms_quadSymPrims) x (dim of atom 9 or 11)

% !! scan Row-wise

function quadSymVec = cVec(quadSymPrims)

[r, c] = size(quadSymPrims);
if ~(r == 1 || c == 1)
    error('Input parameter must be either row or column vector');
end

quadSymVec = [];

for j = 1:c
    for i = 1:r
        
        quadSymPrim = quadSymPrims(i,j);
        
        v = [quadSymPrim.pos' quadSymPrim.r  quadSymPrim.Up1' quadSymPrim.Um1' ...
            quadSymPrim.elongation];
        if (r==1 && c~=1)
            quadSymVec = [quadSymVec; v];
       end
    end
     if (c==1 && r~=1)
        quadSymVec = [quadSymVec; v];
     end  
end

if (r==1 && c ==1)
    quadSymVec = v;
end
