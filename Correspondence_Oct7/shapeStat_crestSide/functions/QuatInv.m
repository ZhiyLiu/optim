function [Qinv] = QuatInv(Q)
% [Qinv] = QuartInv(Q)
%
% Compute the inverse of quarternions.
% Input: - Q = w + x*i + y*j + z*k ; 
% Output:  Qinv = w - x*i - y*j - z*k;

Qinv = [Q(1,:); -Q(2:4,:)];