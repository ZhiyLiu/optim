function [v] = SphereLog(V)
% [v] = SphereLog(V)
%
% Compute the logrithm in S^2.
% Input: - V = (x, y, z) = (cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi));
% Output:  q = theta*(cos(phi), sin(phi)); 

theta = real(acos(V(1,:)));
id1 = find(abs(theta) > 1e-16);

fac = ones(size(theta));
fac(id1) = theta(id1) ./ sin(theta(id1));

v = repmat(fac, 2, 1) .* V(2:3,:);