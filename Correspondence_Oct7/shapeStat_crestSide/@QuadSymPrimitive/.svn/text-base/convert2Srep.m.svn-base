function srepPrims = convert2Srep(qSymPrims)
% quadSimPrims = convert2Srep(srepQuadPrims)
%
% Convert symmetric space primitives to s-rep (Morphormics) primitives.

[r,c] = size(qSymPrims);

srepPrims = [];

for i = 1:r
    aSrepRow = [];
    for j = 1:c
        qSymPrim = qSymPrims(i,j);
        if (qSymPrim.inTangentSpace)
            error('This atom is in the tangent space, not in the manifold');
        end

        r = qSymPrim.r;
        r3 = r * qSymPrim.elongation;
        r = [r r r3];

        v1 = qSymPrim.Up1;
        v2 = qSymPrim.Um1;
        v3 = bisector(v1, v2);
        U = [v1 v2 v3];

        srepPrim = SrepQuadPrimitive(qSymPrim.pos, r, U, false);

        aSrepRow = [aSrepRow , srepPrim];
    end
    srepPrims = [srepPrims;  aSrepRow];
end

function vb = bisector(v1, v2)

vb = v1 + v2;
nb = norm(vb);
if nb == 0
    vb = cross(v1, [0; 1; 0]);
    nb = norm(vb);
    if nb == 0 % In case v1 is collinear with (0, 1, 0)
        vb = [1; 0; 0];
        return;
    end
end

vb = vb / nb;
return
