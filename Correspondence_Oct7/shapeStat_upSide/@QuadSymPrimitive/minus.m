function diffQSymPrim = minus(qSP2, qSP1)

if (qSP2.inTangentSpace)
    if (qSP1.inTangentSpace)
        diffQSymPrim = QuadSymPrimitive(qSP2.pos - qSP1.pos, ...
            qSP2.r - qSP1.r, qSP2.elongation - qSP1.elongation, ...
            qSP2.Up1 - qSP1.Up1, qSP2.Um1-qSP1.Um1, true);
    else
        error('Both operands must be in the same space.');
    end

else
    if (qSP1.inTangentSpace)
        error('Both operands must be in the same space.');
    else
        diffQSymPrim = QuadSymPrimitive(qSP2.pos - qSP1.pos, ...
            qSP2.r / qSP1.r, qSP2.elongation / qSP1.elongation, ...
            symAction(qSP1.Up1, qSP2.Up1), symAction(qSP1.Um1, qSP2.Um1), false);
    end
end



