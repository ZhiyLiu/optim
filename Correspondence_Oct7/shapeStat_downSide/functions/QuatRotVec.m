function [vq] = QuatRotVec(Q, v)
%
% Rotate vectors v by quaternions Q
% Input: - Q = (w, x, y, z) representing w + x*i + y*j + z*k; 
%        - v = (vx, vy, vz);
% Output:  vq = QvQ^{-1}

[nDim1, nq] = size(Q);
[nDim2, nv] = size(v);

if (nq == 1)
    Q = Q * ones(1, nv);
    nq = nv;
elseif (nv == 1)
    v = v * ones(1, nq);
    nv = nq;    
elseif (nq ~= nv)
    disp('Error in QuatRotVec(Q, v): Quats and vecs differ in number!');
    return;
end

v = [zeros(1, nv); v];

Q1 = [Q(1,:).*v(1,:) - Q(2,:).*v(2,:) - Q(3,:).*v(3,:) - Q(4,:).*v(4,:); ...
        Q(1,:).*v(2,:) + Q(2,:).*v(1,:) + Q(3,:).*v(4,:) - Q(4,:).*v(3,:); ...
        Q(1,:).*v(3,:) - Q(2,:).*v(4,:) + Q(3,:).*v(1,:) + Q(4,:).*v(2,:); ...
        Q(1,:).*v(4,:) + Q(2,:).*v(3,:) - Q(3,:).*v(2,:) + Q(4,:).*v(1,:)];
Q2 = QuatInv(Q);
vq = [Q1(1,:).*Q2(1,:) - Q1(2,:).*Q2(2,:) - Q1(3,:).*Q2(3,:) - Q1(4,:).*Q2(4,:); ...
        Q1(1,:).*Q2(2,:) + Q1(2,:).*Q2(1,:) + Q1(3,:).*Q2(4,:) - Q1(4,:).*Q2(3,:); ...
        Q1(1,:).*Q2(3,:) - Q1(2,:).*Q2(4,:) + Q1(3,:).*Q2(1,:) + Q1(4,:).*Q2(2,:); ...
        Q1(1,:).*Q2(4,:) + Q1(2,:).*Q2(3,:) - Q1(3,:).*Q2(2,:) + Q1(4,:).*Q2(1,:)];

vq = vq(2:4, :);
