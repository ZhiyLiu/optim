%==========================================================================
% function qSymPrim = Exp(L)
%
%   Compute the exp map of a primitive of SYMM representation.  
%   L is a matrix whose columns are logs of atoms.  
%   The columns of M are the corresponding exp map.
%--------------------------------------------------------------------------

function [prims] = Exp(tangentPrims)

if( sum([tangentPrims.inTangentSpace]) ~= numel(tangentPrims) )
	tangentPrims
	error( 'symmetric-space Quad Primitive is not in tangent space.' );
end

prims	= tangentPrims;

for i = 1:numel(tangentPrims)
	prim	    = tangentPrims(i);
	prim.r		= exp(prim.r);
	prim.elongation	= exp(prim.elongation);
	prim.Up1	= SphereExp(prim.Up1);
    prim.Um1    = SphereExp(prim.Um1);
	prim.inTangentSpace	= false;
	prims(i)	= prim;
end

return
