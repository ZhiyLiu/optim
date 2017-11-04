function [tubeSymVec] = cVec(prims);
%
% [tubeSymVec] = cVec(prim)
%
% returns a column vector containing the elements of the symmetric space
% tube primitive 'prim' in the order pos, r, U0, elongation
%
% prim can also be a column or a row vector of primitives,
% in which case, the return value is sample size x dim of atom.
%

[r, c]	= size(prims);
if ~(r == 1 | c == 1)
	error('Input parameter prims must be either a row or a column vector');
	return
end

tubeSymVec	= [];

for i = 1:r
	for j = 1:c
		p	= prims(i,j);

		% TODO: Should this be fixed to include dr or not?
		v = [p.pos' p.r p.U0' p.hca p.elongation];
		tubeSymVec	= [tubeSymVec; v];
	end
end

return;

