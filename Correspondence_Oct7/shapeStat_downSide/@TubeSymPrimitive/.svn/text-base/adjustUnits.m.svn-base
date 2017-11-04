function [adjusted] = adjustUnits(prim, r, bCommensurate, varargin)
%
% [adjusted] = adjustUnits(prim, r, bCommensurate)
%
% Adjusts the units of the primitive prim by multiplying (bCommensurate = 1),
% or dividing (bCommensurate = 0) by the factor 'r'.
%

adjusted	= TubeSymPrimitive();
if( prim.inTangentSpace == 1)
	if ( bCommensurate == 0)
		adjusted.pos		= prim.pos;
		adjusted.r			= prim.r / r;
		adjusted.U0			= prim.U0 / r;
		adjusted.hca		= prim.hca / r;
		adjusted.elongation	= prim.elongation / r;
		adjusted.dr			= prim.dr / r;
	else
		adjusted.pos		= prim.pos;
		adjusted.r			= r * prim.r;
		adjusted.U0			= r * prim.U0;
		adjusted.hca		= r * prim.hca;
		adjusted.elongation	= r * prim.elongation;
		adjusted.dr			= r * prim.dr;
	end
	adjusted.inTangentSpace	= prim.inTangentSpace;
	adjusted.baseAtom		= prim.baseAtom;
else
	prim
	error('Primitive prim should be in tangent space.');
end
return;
