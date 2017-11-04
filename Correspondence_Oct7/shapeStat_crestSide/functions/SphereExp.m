%=================================================================================
% Function M = SphereExp(L)
%
%   Compute the exp map on S^2.  L is a matrix whose columns are logs in S^2.
%     The columns of M are the corresponding exp map.
%---------------------------------------------------------------------------------

function [M] = SphereExp(L)

[d1, d2, d3] = size(L);

if (d1 ~= 2)
    uiwait(errordlg('Input dimension in SphereExp is probably wrong!'));
end

% L = (cos(phi), sin(phi)) * theta
L = reshape(L, d1, d2*d3);
theta = sqrt(sum(L.^2));
id1 = find(abs(theta) > 1e-16);

fac = ones(size(theta));
fac(id1) = sin(theta(id1)) ./ theta(id1);

M = [cos(theta); repmat(fac,d1,1).*L];
M = squeeze(reshape(M, d1+1,d2,d3));
