% Use to put data to scale by World Extent
%
% sWE: scale 

function qSymPrim = toScale(qSymPrim, sWE)

qSymPrim = QuadSymPrimitive(sWE*qSymPrim.pos, sWE*qSymPrim.r, ...
                          qSymPrim.elongation, qSymPrim.Up1, qSymPrim.Um1);

