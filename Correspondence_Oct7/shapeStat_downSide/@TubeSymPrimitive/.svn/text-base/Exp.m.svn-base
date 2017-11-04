function [prims] = Exp(tangentPrims);
%
% [prim] = Exp(tangent)
%
% Compute the exponential map of the tube primitive in tangent.
%
% tangent can be an array of tube primitives.
%

if( sum([tangentPrims.inTangentSpace]) ~= prod(size(tangentPrims)) )
	tangentPrims
	error( 'symmetric-space Tube Primitive is not in tangent space.' );
end

prims	= tangentPrims;

for i = 1:prod(size(tangentPrims))
	% FIXME: verify RP1 stats.
	prim	= tangentPrims(i);
	%prim.pos		= tangent.pos;	%already copied
	prim.r			= exp(prim.r);
	prim.elongation	= exp(prim.elongation);
	prim.hca		= atan(prim.hca) + pi/2.0;
	prim.U0			= SphereExp(prim.U0);
	prim.inTangentSpace	= false;
	prims(i)	= prim;
end

return
