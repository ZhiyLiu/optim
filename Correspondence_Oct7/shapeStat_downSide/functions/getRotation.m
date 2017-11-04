function [Q] = getRotation(S)
% [Q] = getRotation(S)
%
% Compute the rotation that takes a point on 
%  S^2 to the point (1,0,0).

% The rotation that takes s = (x, y, z) on S^2 to (1, 0, 0) is given by 
% {|(0, z, -y)|, theta}, where cos(theta) = x.  

theta2 = real(acos(S(1,:)))/2;

%sn = sqrt(sum(S(2:3, :).^2));
id1 = find(theta2 > 1e-16);
%fac = ones(size(theta2));
%fac(id1) = sin(theta2(id1)) ./ sqrt(sum(S(2:3,id1).^2));

fac = ones(size(theta2)) ./ cos(theta2) / 2;
Q = [cos(theta2); zeros(size(theta2)); S(3,:).*fac; -S(2,:).*fac];

