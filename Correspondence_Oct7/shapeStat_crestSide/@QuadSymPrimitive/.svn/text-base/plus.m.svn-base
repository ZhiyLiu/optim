function [qSP] = plus(qSP1, deltaQSP2)

qSP = qSP1;

for i = 1:numel(qSP1)
    pos = qSP1(i).pos + deltaQSP2(i).pos;
    r = qSP1(i).r * deltaQSP2(i).r;
    Up1 = QuatRotVec( QuatInv(getRotation(qSP1(i).Up1)), deltaQSP2(i).Up1);
    Um1 = QuatRotVec( QuatInv(getRotation(qSP1(i).Um1)), deltaQSP2(i).Um1);
    elongation = qSP1(i).elongation * deltaQSP2(i).elongation;

    qSP(i) = QuadSymPrimitive(pos,r,elongation, Up1, Um1, false);
end



