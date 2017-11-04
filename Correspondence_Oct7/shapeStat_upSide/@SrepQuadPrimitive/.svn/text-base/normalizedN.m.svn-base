function normN = normalizedN(qSrepPrim)

if (qSrepPrim.inTangentSpace)
    error('This atom is already in Tangent Space.');
else
    v1 = qSrepPrim.U(:,1);
    v2 = qSrepPrim.U(:,2);
    n1 = norm(v1);
    n2 = norm(v2);
    if n1 == 0 || n2 == 0
        normN = zeros(3, 1);
    end
    normN = v1/n1 - v2/n2;
    normN = normN / norm(normN);
end

return;
