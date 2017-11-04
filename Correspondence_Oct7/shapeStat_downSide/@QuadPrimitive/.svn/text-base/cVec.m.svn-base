% function v = cVec(qPrim) returns a column vector of quadPrimitive 'qPrim'
%   

% qPrims is ONE or a ROW or COLUMN vector of QuadPrimitives
%
% Output quadVec : (probably sample size) x  (dim of atom 9 or 10)

function quadVec = cVec(qPrims)

[r, c] = size(qPrims);
if ~(r == 1 || c == 1)
    error('Input parameter must be either row or column vector');
end

quadVec = [];

for i = 1:r
    for j = 1:c
        
       qPrim = qPrims(i,j);
        
       v = [qPrim.pos' qPrim.r qPrim.q' qPrim.theta qPrim.elongation];
       % position, radius, quaternion, elongation
       if (r==1 && c~=1)
            quadVec = [quadVec; v];
       end
    end
    if (c==1 && r~=1)
        quadVec = [quadVec; v];
    end
end

if (r==1 && c ==1)
    quadVec = v;
end
