function [val] = get(prim, propName);
%
% [val] = get(prim, propName)
%
% Returns a property of the symmetric-space tube primitive prim by the name
% of propName.
%
% Possible property names are
% x,y,z, pos, r, U0, elongation, theta (hca), isTangentSpace
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
	case 'U0'
		val = [prim.U0];
	case 'elongation'
		val = [prim.elongation];
	case 'theta'
		val = [prim.hca];
	case 'hca'
		val = [prim.hca];
	case 'baseAtom'
		val = [prim.baseAtom];
	case 'dr'
		val	= [prim.dr];
	case 'inTangentSpace'
		val = [prim.inTangentSpace];
	otherwise
		error([propName ' is not a valid TubeSymPrimitive property'])
end

return;
