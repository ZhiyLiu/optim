function d = squaredNorm(SrepPrim)

% d = squaredNorm(SrepPrim)
% 
% Treats every coefficient as an entry in a big vector and computes the
% squared norm.  For use in the tangent space (?)  Gives
% a different result than squaredNorm(QuadSymPrim) because the coefficients
% are different.

v = [SrepPrim.pos' SrepPrim.r reshape(SrepPrim.U, 1, [])];
d = sum(v .* v);
