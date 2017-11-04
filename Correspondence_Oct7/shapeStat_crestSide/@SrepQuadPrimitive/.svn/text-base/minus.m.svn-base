function diffSrPrim = minus(srp2, srp1)

inTangentSpace = srp2.inTangentSpace;

% Error if truth of "inTangentSpace" differs
if ~srp2.inTangentSpace ~= ~srp1.inTangentSpace
    error('Both operands must be in the same space.');
end

if inTangentSpace
    disp('hmm')
    diffSrPrim = SrepQuadPrimitive(srp2.pos - srp1.pos, ...
          srp2.r - srp1.r, srp2.U - srp1.U, true);
else
    diffSrPrim = SrepQuadPrimitive(srp2.pos - srp1.pos, ...
          srp2.r ./ srp1.r, ...
          [ symAction(srp1.U(:, 1), srp2.U(:, 1)) ...
              symAction(srp1.U(:, 2), srp2.U(:, 2)) ...
              symAction(srp1.U(:, 3), srp2.U(:, 3))], false);
end



