function [tubePrims] = convert2Lie(symPrims);
%
% [prim] = TubePrimitive(symPrims);
%
% Converts the tube primitive symPrims from symmetric space representation to
% lie group representation.
%
% symPrims can also be an array of atoms.
%

[r,c]	= size(symPrims);
tubePrims	= [];

for i = 1:r
	aTubeRow	= [];
	for j = 1:c
		sym		= symPrims(i,j);
		
		% FIXME: Rohit Need to figure out a more stable way to find out Up1 and Um1.
		% It doesn't particularly matter what we do here, as long as
		% alignTube is called on the entire tube in the end.
		perp	= cross( sym.U0, [ 1; 0; 0.01 ]/sqrt(1+0.01*0.01) );
		if( sum(perp.*perp) <= 1e-8 )
			perp	= cross( sym.U0, [0; 1; 0 ]);
		end
		perp	= perp/sqrt(sum(perp.*perp));
		cost	= repmat(cos(sym.hca), [ 3 1 ] );
		sint	= repmat(sin(sym.hca), [ 3 1 ] );
		Up1		= sym.U0.*cost + perp.*sint;
		Um1		= sym.U0.*cost - perp.*sint;

		q	= quatFromFrame( sym.U0, perp );
        
		prim	= TubePrimitive(sym.pos, sym.r, sym.elongation, q, sym.hca, false, sym.baseAtom, sym.dr );

		%
		% End of conversion
		%
		aTubeRow	= [ aTubeRow, prim ];
	end

	tubePrims	= [tubePrims; aTubeRow ];
end

return;
