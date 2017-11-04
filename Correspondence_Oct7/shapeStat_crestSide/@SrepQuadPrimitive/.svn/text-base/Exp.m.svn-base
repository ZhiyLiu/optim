%==========================================================================
% function qSymPrim = Exp(L)
%
%   Compute the exp map of a primitive in the s-rep (Mx) representation.  
%   L is a structure array whose entries are logs of quad primitives (atoms)
%   The entries of M are the corresponding quad primitives on the manifold.
%--------------------------------------------------------------------------

function [prims] = Exp(tangentPrims)

if( sum([tangentPrims.inTangentSpace]) ~= numel(tangentPrims) )
    tangentPrims
    error( 'symmetric-space Quad Primitive is not in tangent space.' );
end

prims	= tangentPrims;

for n = 1:numel(tangentPrims)
    prim = tangentPrims(n);
    prim.r = exp(prim.r);
    prim.U = SphereExp(prim.U);
    prim.inTangentSpace = false;
    prims(n) = prim;
end

return
