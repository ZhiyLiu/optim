function SrepQuadPrims = convert2Srep(qPrims)

% srepQuadPrims = convert2Srep(qPrims)
%  qPrims: a structure array of quad primitives (or just a single one)

[r,c] = size(qPrims);

SrepQuadPrims = [];

for i = 1:r
    aSrepQuadRow = [];
    for j = 1:c
        qPrim = qPrims(i,j);
        if (qPrim.inTangentSpace)
            error('This atom is in the tangent space, not the manifold');
        end

        cos_t = cos(qPrim.theta);
        sin_t = sin(qPrim.theta);

        r = [qPrim.r qPrim.r qPrim.r * qPrim.elongation];

        v1 = QuatRotVec(qPrim.q, [cos_t; sin_t; 0]);
        v2 = QuatRotVec(qPrim.q, [cos_t; -sin_t; 0]);
        v3 = QuatRotVec(qPrim.q, [1; 0; 0]);
        U = [v1 v2 v3];

        SrepQuadPrim = SrepQuadPrimitive(qPrim.pos, r, U, false);

        aSrepQuadRow = [aSrepQuadRow , SrepQuadPrim];
    end
      SrepQuadPrims = [SrepQuadPrims;  aSrepQuadRow];
end
