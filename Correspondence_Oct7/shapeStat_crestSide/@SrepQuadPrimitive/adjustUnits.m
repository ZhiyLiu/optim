function qPrim = adjustUnits(qPrim, mr, bCommensurate, weightMode)

% bCommensurate -> 1 : multiplication, bCommensurate -> 0 : division 
% mr: mean radius 

% Make sure it's in the tangent space.
if( qPrim.inTangentSpace ~= 1 )
    qPrim
    error('QuadPrimitive prim should be in tangent space.');
end

if (nargin > 3 && weightMode) 
    w = mr/2;
else
    w = mr;
end

if (bCommensurate)
    qPrim = SrepQuadPrimitive(qPrim.pos, qPrim.r*mr, qPrim.U*w, true);
else
    qPrim = SrepQuadPrimitive(qPrim.pos, qPrim.r/mr, qPrim.U/w, true);
end

return;









