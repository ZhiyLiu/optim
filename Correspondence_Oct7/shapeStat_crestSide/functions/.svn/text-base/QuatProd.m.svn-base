% FUNCTION [Q] = QuatProd(Q1, Q2)
%
% Compute quaternion products
% Input: - Q1 = (w1, x1, y1, z1) representing w1 + x1*i + y1*j + z1*k; 
%        - Q2 = (w2, x2, y2, z2) representing w2 + x2*i + y2*j + z2*k; 
% Output:  Q = Q1.*Q2

function [Q] = QuatProd(Q1, Q2)

[nDim1, nq1] = size(Q1);
[nDim2, nq2] = size(Q2);

if (nq1 == 1)
    Q1 = Q1 * ones(1, nq2);
    nq1 = nq2;
elseif (nq2 == 1)
    Q2 = Q2 * ones(1, nq1);
    nq2 = nq1;    
elseif (nq1 ~= nq2 | nDim1 ~= nDim2)
    disp('Error in QuatProd(Q1, Q2): Q1 and Q2 must have the same size!');
    return;
end

Q = [Q1(1,:).*Q2(1,:) - Q1(2,:).*Q2(2,:) - Q1(3,:).*Q2(3,:) - Q1(4,:).*Q2(4,:); ...
        Q1(1,:).*Q2(2,:) + Q1(2,:).*Q2(1,:) + Q1(3,:).*Q2(4,:) - Q1(4,:).*Q2(3,:); ...
        Q1(1,:).*Q2(3,:) - Q1(2,:).*Q2(4,:) + Q1(3,:).*Q2(1,:) + Q1(4,:).*Q2(2,:); ...
        Q1(1,:).*Q2(4,:) + Q1(2,:).*Q2(3,:) - Q1(3,:).*Q2(2,:) + Q1(4,:).*Q2(1,:)];

% make sure w>=0
cols = find(Q(1,:) < 0);
Q(:, cols) = -Q(:, cols);