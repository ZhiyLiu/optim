%==========================================================================
% function qPrim = Exp(logQPrim)
%
%   Compute the exp map of a primitive of LIE gp representation.  L is a matrix whose columns are logs of
%   atoms.  The columns of M are the corresponding exp map.
%--------------------------------------------------------------------------

function qPrim = Exp(logQPrim)

if (~logQPrim.inTangentSpace)
    error('This atom is in manifold')
    return;
end

qPrim = QuadPrimitive(logQPrim.pos, exp(logQPrim.r), ...
    exp(logQPrim.elongation), QuatExp(logQPrim.q), logQPrim.theta, false);
