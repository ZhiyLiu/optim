function [q] = QuatLog(Q)
% [q] = QuatLog(Q)
%
% Compute the exponential for quarternions.
% Input: - Q = w + x*i + y*j + z*k = cos(theta/2) + sin(theta/2)*v;
% Output:  q = theta*(vx, vy, vz); 

theta = real(acos(Q(1,:))) * 2;

id1 = find(abs(theta)>1e-8);
fac = ones(size(theta));
fac(id1) = theta(id1)./sin(theta(id1)/2);

q = [Q(2:4,:) .* repmat(fac, 3,1)];