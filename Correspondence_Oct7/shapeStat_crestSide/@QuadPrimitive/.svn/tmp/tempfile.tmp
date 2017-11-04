function [] = display(p)
%
% function [] = display(p)
%
% Displays the primitive value in Lie group space.
%

[r, c] = size(p);
for i = 1:r
    for j = 1:c
        quadPrim = p(i,j);
		disp([ 'QuadPrimitive (Lie group space) [' num2str(i) ', ' num2str(j) ']' ]);
        disp([ '  position         : (' num2str(quadPrim.pos(1)) ', ' num2str(quadPrim.pos(2)) ', ' num2str(quadPrim.pos(3)) ')' ]);
        disp([ '  r                : ' num2str(quadPrim.r) ]);
        disp([ '  elongation       : ' num2str(quadPrim.elongation) ]);
        disp([ '  quat             : (' num2str(quadPrim.q(1)) ', ' num2str(quadPrim.q(2)) ', ' num2str(quadPrim.q(3)) ', ' num2str(quadPrim.q(4)) ')' ]);
        disp([ '  theta(degrees)   : ' num2str(quadPrim.theta * 180/pi) ]);
        disp([ '  in tangent space : ' int2str(quadPrim.inTangentSpace) ]);
    end
end
return;
