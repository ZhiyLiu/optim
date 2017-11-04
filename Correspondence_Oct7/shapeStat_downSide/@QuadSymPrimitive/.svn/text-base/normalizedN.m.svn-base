function normN = normalizedN(qSymPrim)

if (qSymPrim.inTangentSpace)
    error('This atom is already in Tangent Space.');
else
    normN = qSymPrim.Up1 - qSymPrim.Um1;
    normN = normN / sqrt(sum(normN.*normN));
end

return;