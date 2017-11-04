% ================================================
% FUNCTION xformedQP = applyMrepXform(qP, T)
%
%   Apply a global transformation T to qP
%   T = [Q; R] where Q is a rotation (quaternion)
%       and R is a scaling factor
% ------------------------------------------------

function xformedQP = applyMrepXform(qP, T)


R = T(5);
Q = T(1:4);

pos = R*QuatRotVec(Q, qP.pos);
r = R*qP.r;
q = QuatProd(Q, qP.q);

xformedQP = QuadPrimitive(pos, r, qP.elongation, q, qP.theta, false);



