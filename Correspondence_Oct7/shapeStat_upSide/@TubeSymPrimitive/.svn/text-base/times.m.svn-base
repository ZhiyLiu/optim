function [prims] = times( prims, scalar );
%
% [t] = prim * scalar;
% 
% Simply multiplies all the values in the prim with the scalar scalar.
% prim can be an array.
%

for i = 1:prod(size(prims))
	prim	= prims(i);
	prim.pos	= prim.pos * scalar(i);
	prim.r		= prim.r * scalar(i);
% TODO: Temporary FIX
%	prim.dr		= prim.r * scalar(i);
	prim.dr		= prim.r;
	prim.elongation	= prim.elongation * scalar(i);
	prim.hca	= prim.hca * scalar(i);
	prim.U0		= prim.U0 * scalar(i);
	prims(i)	= prim;
end

return;

