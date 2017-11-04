function srp = plus(srp1, deltaSrp)

srp = srp1;

for i = 1:numel(srp1)
    pos = srp1(i).pos + deltaSrp(i).pos;
    r = srp1(i).r .* deltaSrp(i).r;
    for n = 1:3
        U1 = srp1(i).U(:,n);
        deltaU = deltaSrp(i).U(:,n);
        U(:,n) = QuatRotVec( QuatInv(getRotation(U1)), deltaU);
    end
    srp(i) = SrepQuadPrimitive(pos, r, U, false);
end
