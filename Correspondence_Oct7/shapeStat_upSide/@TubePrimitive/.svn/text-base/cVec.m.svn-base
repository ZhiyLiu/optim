% function v = cVec(tPrim) returns a column vector of tubePrimitive 'tPrim'
%   

% tPrims is ONE or a ROW or COLUMN vector of TubePrimitives
%
% Output tubeVec : (probably sample size) x  (dim of atom 9 or 10)

function tubeVec = cVec(tPrims)

[r, c] = size(tPrims);
if ~(r == 1 | c == 1)
    error('Input parameter must be either row or column vector');
    return 
end

tubeVec = [];

for i = 1:r
    for j = 1:c
        
       tPrim = tPrims(i,j);
        
	   % TODO: To include or not to include spokes.
       v = [tPrim.pos' tPrim.r tPrim.q' tPrim.theta tPrim.elongation];
       % position, radius, quaternion, elongation
       if (r==1 && c~=1)
            tubeVec = [tubeVec; v];
       end
    end
    if (c==1 && r~=1)
        tubeVec = [tubeVec; v];
    end
end

if (r==1 && c ==1)
    tubeVec = v;
end
