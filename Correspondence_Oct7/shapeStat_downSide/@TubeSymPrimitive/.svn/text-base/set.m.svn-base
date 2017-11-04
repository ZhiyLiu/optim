function [prim] = set( prim, varargin );
%
% [prim] = set(prim, [propName, value]* )
%
% Sets a property of the symmetric-space tube primitive prim by the name of
% propName to value.
%
% Possible property names are
% x,y,z, pos, r, U0, elongation, theta (hca), isTangentSpace
%

for i=1:2:length(varargin)
	val	= varargin{i+1};
	switch (varargin{i})
		case 'x'
			prim.pos(1)	= val;
		case 'y'
			prim.pos(2)	= val;
		case 'z'
			prim.pos(3)	= val;
		case 'pos'
			prim.pos	= val;
		case 'r'
			prim.r		= val;
		case 'U0'
			prim.U0		= val;
		case 'elongation'
			prim.elongation		= val;
		case 'theta'
			prim.hca	= val;
		case 'hca'
			prim.hca	= val;
		case 'baseAtom'
			prim.baseAtom	= val;
		case 'dr'
			prim.dr		= val;
		case 'inTangentSpace'
			prim.inTangentSpace	= val;
		otherwise
			error([varargin{i} ' is not a valid TubeSymPrimitive property'])
	end
end

return;
