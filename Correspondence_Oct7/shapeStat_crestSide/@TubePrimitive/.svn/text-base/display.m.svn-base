function [] = display(prims)
%
% function [] = display(prims)
%
% Displays the primitive value in Lie group space.
%

[r, c] = size(prims);
for i = 1:r
    for j = 1:c
        p = prims(i,j);
		disp([ 'TubePrimitive (Lie group space) [' num2str(i) ', ' num2str(j) ']' ]);
		disp([ '  position         : (' num2str(p.pos(1)) ', ' num2str(p.pos(2)) ', ' num2str(p.pos(3)) ')' ]);
		disp([ '  r                : ' num2str(p.r) ]);
		disp([ '  elongation       : ' num2str(p.elongation) ]);
		disp([ '  quat             : (' num2str(p.q(1)) ', ' num2str(p.q(2)) ', ' num2str(p.q(3)) ', ' num2str(p.q(4)) ')' ]);
		disp([ '  theta(degrees)   : ' num2str(p.theta * 180/pi) ]);
		disp([ '  in tangent space : ' int2str(p.inTangentSpace) ]);
		disp([ '  base Atom        : ' int2str(p.baseAtom) ]);
		disp([ '  dr               : (' num2str(p.dr') ')' ]);
	end
end
return;
