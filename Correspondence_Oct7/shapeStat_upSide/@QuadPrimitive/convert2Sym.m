function quadSymPrims = convert2Sym(qPrims)

[r,c] = size(qPrims);

quadSymPrims = [];

for i = 1:r
    aQuadSymRow = [];
    for j = 1:c
        qPrim = qPrims(i,j);
        if (qPrim.inTangentSpace)
            error('This atom is in tangent space, not in manifold');
        end

        cos_t = cos(qPrim.theta);
        sin_t = sin(qPrim.theta);

        v1 = QuatRotVec(qPrim.q, [cos_t; sin_t; 0]);
        v2 = QuatRotVec(qPrim.q, [cos_t; -sin_t; 0]);

        quadSymPrim = QuadSymPrimitive(qPrim.pos, qPrim.r, qPrim.elongation, ...
            v1, v2, false);

        aQuadSymRow = [aQuadSymRow , quadSymPrim];
    end
      quadSymPrims = [quadSymPrims;  aQuadSymRow];
end
