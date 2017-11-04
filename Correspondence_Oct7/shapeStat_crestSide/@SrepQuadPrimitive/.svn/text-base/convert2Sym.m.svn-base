function quadSymPrims = convert2Sym(qSrepPrims)
% quadSimPrims = convert2Sym(srepQuadPrims)
%
% Convert Mx-format (s-rep) primitives to symmetric representation
% primitives in the UNC format.  This is a lossy conversion since the
% opposing spokes may have different radii in the s-rep representation 
% but are forced to have the same radius in the UNC format.  The geometric 
% mean is used if they are different in the input.

[r,c] = size(qSrepPrims);

quadSymPrims = [];

for i = 1:r
    aQuadSymRow = [];
    for j = 1:c
        qSrepPrim = qSrepPrims(i,j);
        if (qSrepPrim.inTangentSpace)
            error('This atom is in tangent space, not in manifold');
        end

        inr = qSrepPrim.r;
        r = sqrt(inr(1) * inr(2));
        elong = inr(3) / r;
        v1 = qSrepPrim.U(:,1);
        v2 = qSrepPrim.U(:,2);

        quadSymPrim = QuadSymPrimitive(qSrepPrim.pos, r, elong, v1, v2, false);

        aQuadSymRow = [aQuadSymRow , quadSymPrim];
    end
    quadSymPrims = [quadSymPrims;  aQuadSymRow];
end
