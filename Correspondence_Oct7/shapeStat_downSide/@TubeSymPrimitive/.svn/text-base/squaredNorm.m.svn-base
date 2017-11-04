function d = squaredNorm(symPrim)

%
% TODO: Include dr here or not?
%

d = sum((symPrim.pos).*(symPrim.pos)) + symPrim.r * symPrim.r + ...
    sum((symPrim.U0) .* (symPrim.U0)) + symPrim.hca.*symPrim.hca + ...
    symPrim.elongation * symPrim.elongation ;
