% ================================================
% FUNCTION xformedQP = applyMrepXform(prim, T)
%
%   Apply a global transformation T to prim
%   T = [Q; R] where Q is a rotation (quaternion)
%       and R is a scaling factor
% ------------------------------------------------

function xformedPrim = applyMrepXform(prim, T)


R = T(5);
Q = T(1:4);

pos = R*QuatRotVec(Q, prim.pos);
r = R*prim.r;
q = QuatProd(Q, prim.q);

xformedPrim = TubePrimitive(pos, r, prim.elongation, q, prim.theta, false, prim.baseAtom, prim.dr);



