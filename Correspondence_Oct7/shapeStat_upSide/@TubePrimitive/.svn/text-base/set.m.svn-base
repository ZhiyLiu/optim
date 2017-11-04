function [prim] = set( prim, varargin );
%
% [prim] = set(prim, [propName, value]* )
%
% Sets a property of the tube primitive prim by the name of propName to value.
%
% Possible property names are
% x,y,z, pos, r, q, elongation, theta, inTangentSpace
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
		case 'q'
			prim.q		= val;
		case 'elongation'
			prim.elongation	= val;
		case 'theta'
			prim.theta		= val;
		case 'baseAtom'
			prim.baseAtom	= val;
		case 'inTangentSpace'
			prim.inTangentSpace	= val;
		case 'dr'
			prim.dr		= val;
		otherwise
			error([varargin{i} ' is not a valid TubePrimitive property'])
	end
end

return;
