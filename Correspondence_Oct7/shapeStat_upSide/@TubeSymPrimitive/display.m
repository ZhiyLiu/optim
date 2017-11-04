function [] = display(sP)
%
% function [] = display(sP)
%
% Displays the primitive value in symmetric space.
%

[r,c]	= size(sP);

for i = 1:r
	for j = 1:c
		p	= sP(i,j);
		disp([ 'TubePrimitive (symmetric space) [' num2str(i) ', ' num2str(j) ']']);
		disp([ '  position         : (' num2str(p.pos(1)) ', ' num2str(p.pos(2)) ', ' num2str(p.pos(3)) ')' ]);
		disp([ '  r                : ' num2str(p.r) ]);
		disp([ '  elongation       : ' num2str(p.elongation) ]);
		if( p.inTangentSpace )
			disp([ '  U0               : (' num2str(p.U0(1)) ', ' num2str(p.U0(2)) ')' ]);
		else
			disp([ '  U0               : (' num2str(p.U0(1)) ', ' num2str(p.U0(2)) ', ' num2str(p.U0(3)) ')' ]);
		end
		disp([ '  theta(degrees)   : ' num2str(p.hca * 180/pi) ]);
		disp([ '  base atom        : ' int2str(p.baseAtom) ]);
		disp([ '  in tangent space : ' int2str(p.inTangentSpace) ]);
		disp([ '  spoke deviations : ' num2str(p.dr') ]);
	end
end

return;
