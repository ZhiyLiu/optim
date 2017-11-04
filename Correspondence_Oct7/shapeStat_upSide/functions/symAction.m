% ==================================================================================
% Function [s] = symAction(s1, s2)
%
%   Computes the action of s1 on s2.  
%   Inputs:     s1, s2 are two matrices, the columns of which are unit 3-vectors.
% ----------------------------------------------------------------------------------

function [s] = symAction(s1, s2)

[nDims1, ns1, na1] = size(s1);
[nDims2, ns2, na2] = size(s2);
ns = ns1; na = na1;

natoms1 = ns1*na1;
natoms2 = ns2*na2;
s1 = reshape(s1, nDims1, natoms1);
s2 = reshape(s2, nDims2, natoms2);

if (natoms1 == 1)
    s1 = s1 * ones(1, natoms2);
    natoms1 = natoms2;
    ns = ns2; na = na2;
elseif (natoms2 == 1)
    s2 = s2 * ones(1, natoms1);
    natoms2 = natoms1;    
elseif (natoms1 ~= natoms2 || nDims1 ~= nDims2)
    disp('Error in symAction(s1, s2): s1 and s2 must have the same size!');
    return;
end

s = QuatRotVec( getRotation(s1), s2 );
% s2Log = SphereLog(s2);
% s1Log = SphereLog(s1);
% sLog = s2Log - s1Log;
%s = SphereExp(sLog);

s = squeeze(reshape(s, nDims1, ns, na));
