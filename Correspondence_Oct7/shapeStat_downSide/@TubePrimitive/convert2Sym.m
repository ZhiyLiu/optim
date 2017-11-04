function [tubeSymPrims] = convert2Sym(prims)
%
% [S] = convert2Sym(prims)
%
% Converts the tube atom in prim to its symmetric representation and
% returns it in S.
% prims can also be an array of primitives.
%


[r,c]	= size(prims);

tubeSymPrims	= [];

for i = 1:r
	aTubeSymRow	= [];
	for j = 1:c
		prim	= prims(i,j);
		% rotate vector [1,0,0] (the tangent to the medial axis) into
		% position.
		U0	= QuatRotVec( prim.q, [1; 0; 0] );
		S	= TubeSymPrimitive( prim.pos, prim.r, prim.elongation, U0, prim.theta, prim.baseAtom, prim.dr, false );
		aTubeSymRow	= [aTubeSymRow, S];
	end
	tubeSymPrims	= [tubeSymPrims; aTubeSymRow];
end

return;
