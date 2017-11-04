%==========================================================================
% function logQPrim = Log(qPrim)
%
%   Compute the log map of atoms of LIE group representation. 
%
%--------------------------------------------------------------------------

function logQPrim = Log(qPrim)

if ( qPrim.inTangentSpace )
    error('This atom is already in a tangent space');
    return;
end

logQPrim = QuadPrimitive(qPrim.pos, log(qPrim.r), log(qPrim.elongation),...
   QuatLog(qPrim.q), qPrim.theta, true);

