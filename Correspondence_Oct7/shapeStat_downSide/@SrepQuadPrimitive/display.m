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
        pos = quadPrim.pos;
        r = quadPrim.r;
        U = quadPrim.U;
        disp([ 'QuadPrimitive (s-rep space) [' num2str(i) ', ' num2str(j) ']' ]);
        fprintf(1, '  pos              : (%.4g, %.4g, %.4g)\n', ...
              pos(1), pos(2), pos(3));
        
        for k = 1:length(r)
            fprintf(1, '  r(%.4g), U(%.4g)       : %.4g, (%.4g', ...
                  k, k, r(k), U(1, k));
            for n = 2:size(U,1)
                fprintf(1, ', %.4g', U(n, k));
            end
            fprintf(1, ')\n');
        end
        disp([ '  in tangent space : ' int2str(quadPrim.inTangentSpace) ]);
    end
end
return;
