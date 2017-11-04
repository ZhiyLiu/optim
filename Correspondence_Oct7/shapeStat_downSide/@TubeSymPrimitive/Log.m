function [tangentPrims] = Log(prims);
%
% [tangent] = Log(prims)
%
% Compute the tangent map of the tube primitive in prims.
%

if( sum([prims.inTangentSpace]) ~= 0 )
	prims
	error( 'symmetric-space Tube Primitive is not in manifold space.' );
end

% FIXME: verify RP1 stats.
tangentPrims	= prims;

for i = 1:prod(size(prims))
	tangent		= prims(i);
	%tangent.pos		= prim.pos;	% already copied
	tangent.r			= log(tangent.r);
	tangent.elongation	= log(tangent.elongation);
	tangent.hca			= tan(tangent.hca - pi/2.0);
	tangent.U0			= SphereLog( tangent.U0 );
	tangent.inTangentSpace	= 1;
	tangentPrims(i)		= tangent;
end

return
