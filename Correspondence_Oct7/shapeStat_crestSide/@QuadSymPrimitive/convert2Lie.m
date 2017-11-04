function quadPrims = convert2Lie(qSymPrims)

[r,c] = size(qSymPrims);
quadPrims = [];

for i = 1:r
    
    aQuadRow = [];
    for j = 1:c
    
        qSymPrim = qSymPrims(i,j);
        %==================================================================
        % MAIN ROUTINE for conversion
        %==================================================================
        
        Up1 = qSymPrim.Up1;
        Um1 = qSymPrim.Um1;

        R = reshape(eye(3), [9, 1]);

        %%% BISECTOR
        b = (Up1 + Um1); 
        b_norm = sqrt(sum(b.^2));
        %----------------------------------------------------------
        % (Case 1) Handle when Up1 is almost opposite to Um1 (collinear)
        %----------------------------------------------------------
        if ( b_norm < 1e-8 )
            Y = [0;1;0];
            theta2 = real(acos(dot(Up1, Y)))/2;
            axis = [Up1(3); 0; -Up1(1)];
            axis_norm = sqrt(sum(axis.^2));
            if ( axis_norm < 1e-8) % Up1 == Y axis
                axis = [1;0;0];
            else
                axis = axis ./ sqrt(sum(axis.^2));
            end
            Q = [cos(theta2); sin(theta2)*axis];
            b = QuatRotVec(Q, [1;0;0]);
        else
            b = b ./ b_norm;
        end
        %----------------------------------------------------------
        

        %%% NORMAL (to the medial sheet)
        n = (Up1 - Um1);
        n_norm = sqrt(sum(n.^2));
        %----------------------------------------------------------
        % (Case 2) Handle when Up1 is very close to Um1 (collinear)
        %----------------------------------------------------------
        if ( n_norm < 1e-8 )  % Up1 and Um1 are collinear %1e-15
            Q = getRotation(b);  % rotation b to (1, 0, 0)
            n = QuatRotVec(QuatInv(Q), [0; 1; 0]);
        else
            n = n ./n_norm;
        end
        %---------------------------------------------------------
        
		q	= quatFromFrame( b, n );
        
        theta = real(acos(dot(Up1, b)));
        quadPrim = QuadPrimitive(qSymPrim.pos, qSymPrim.r, qSymPrim.elongation, ...
            q, theta, false);

        %==================================================================
        % END Of MAIN ROUTINE for conversion
        %==================================================================
        
        aQuadRow = [aQuadRow, quadPrim];
    end
    
    quadPrims = [quadPrims ; aQuadRow];
end
