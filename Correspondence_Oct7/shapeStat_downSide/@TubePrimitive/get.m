function [val] = get(prim, propName);
%
% [val] = get(prim, propName)
%
% Returns a property of the tube primitive prim by the name of propName.
% prim can be an array of tube primitives, in which case a row array of
% values is returned.
%
% Possible property names are
% x,y,z, pos, r, q, elongation, theta, inTangentSpace
%

switch propName
	case 'x'
		val = [prim.pos(1)];
	case 'y'
		val = [prim.pos(2)];
	case 'z'
		val = [prim.pos(3)];
	case 'pos'
		val	= [prim.pos];
	case 'r'
		val = [prim.r];
	case 'q'
		val = [prim.q];
	case 'elongation'
		val = [prim.elongation];
	case 'theta'
		val = [prim.theta];
	case 'baseAtom'
		val = [prim.baseAtom];
	case 'inTangentSpace'
		val = [prim.inTangentSpace];
	case 'dr'
		val	= [prim.dr];
	otherwise
		error([propName ' is not a valid TubePrimitive property'])
end
return;
