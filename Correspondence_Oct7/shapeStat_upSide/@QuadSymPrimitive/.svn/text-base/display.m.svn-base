function display(sP)
%
% function [] = display(sP)
%
% Displays the primitive value in symmetric space.
%

[r, c] = size(sP);
for i = 1:r
    for j = 1:c
        sPrim = sP(i,j);
		disp([ 'QuadPrimitive (symmetric space) [' num2str(i) ', ' num2str(j) ']']);
        disp([ '  position         : (' num2str(sPrim.pos(1)) ', ' num2str(sPrim.pos(2)) ', ' num2str(sPrim.pos(3)) ')' ]);
        disp([ '  r                : ' num2str(sPrim.r) ]);
        disp([ '  elongation       : ' num2str(sPrim.elongation) ]);
        if (sPrim.inTangentSpace)
            disp([ '  Up1               : (' num2str(sPrim.Up1(1)) ', ' num2str(sPrim.Up1(2))  ')' ]);
            disp([ '  Um1               : (' num2str(sPrim.Um1(1)) ', ' num2str(sPrim.Um1(2))  ')' ]);
        else
            disp([ '  Up1               : (' num2str(sPrim.Up1(1)) ', ' num2str(sPrim.Up1(2)) ', ' num2str(sPrim.Up1(3)) ')' ]);
            disp([ '  Um1               : (' num2str(sPrim.Um1(1)) ', ' num2str(sPrim.Um1(2)) ', ' num2str(sPrim.Um1(3)) ')' ]);
        end
        disp([ '  in tangent space : ' int2str(sPrim.inTangentSpace) ]);
    end
end

return;
