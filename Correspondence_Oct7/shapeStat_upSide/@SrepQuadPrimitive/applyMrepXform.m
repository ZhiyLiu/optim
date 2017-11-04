% ================================================
% FUNCTION xformedQP = applyMrepXform(sP, T)
%
%   Apply a global transformation T to sP
%   T = [Q; R] where Q is a rotation (quaternion)
%       and R is a scaling factor
%   Alternatively, we may have
%   T = [x; Q; R], where x is a translation
%   NOTE: I have been confused about the order of the components when T is
%   included.  This appears to be the convention used by readM3d. -MSF
% ------------------------------------------------

function xformedQP = applyMrepXform(sP, T)

if numel(T) == 5
    Translation = zeros(3,1);
    Q = T(1:4);
    R = T(5);
elseif numel(T) == 8
    Translation = T(1:3);
    Q = T(4:7);
    R = T(8);
end

pos = R*QuatRotVec(Q, sP.pos) + Translation;
r = R*sP.r;
U = QuatRotVec(Q, sP.U);

xformedQP = SrepQuadPrimitive(pos, r, U, false);



