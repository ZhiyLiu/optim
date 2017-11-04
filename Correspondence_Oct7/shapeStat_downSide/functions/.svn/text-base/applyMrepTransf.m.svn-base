% =========================================================================
% FUNCTION [M] = applyMrepTransf(T, M0)
%
%   INPUT: M0 - cell array of quadPrimitives(lie group representation)
%   OUTPUT: M - cell array of quadPrimitives(lie group representation)
%               after applying a global transformation T to M0
%
%   T = [q; r] where q is a rotation (quaternion)
%       and r is a scaling factor
% -------------------------------------------------------------------------

function [M] = applyMrepTransf(T, M0)

nAtoms = length(M0);

for i = 1:nAtoms
    qP = M0{i};
    xformedQP = applyMrepXform(qP, T);
    M{i} = xformedQP;
end

