function [Q] = QuatExp(q)
% [Q] = QuatLog(q)
%
% Compute the exponential for quarternions.
% Input: - q = theta*(vx, vy, vz); 
% Output:  Q = w + x*i + y*j + z*k = cos(theta/2) + sin(theta/2)*v;

theta = sqrt(sum(q.*q));

th = 1e-4;
id2 = find(theta >= th);

Q = repmat([1;0;0;0], 1, length(theta));
if (~isempty(id2))
    Q(:, id2) = [cos(theta(id2)/2); q(:, id2) .* ...
        repmat(sin(theta(id2)/2)./theta(id2), 3,1)];
end

return