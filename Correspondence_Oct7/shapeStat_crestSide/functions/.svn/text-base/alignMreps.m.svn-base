% ================================================
% FUNCTION [T] = alignMreps(M1, M2, alignMode)
%
%   Align mrep M1 to mrep M2, using M2's radii info
%
%   alignMode: 1 - translation 
%              2 - rotation
%              3 - rotation and scaling
% ------------------------------------------------

function [T1] = alignMreps(M1, M2, tMode)

T1 = [];
if tMode ~= 2 && tMode ~=3
    uiwait(errorbox('Unknown alignment mode in alignMreps.  Abort.'));
    return;
end

Tnow = [1;0;0;0;1];
% rv = rand(5,1);
% rv(1) = rv(1)+50;
% rv(1:4) = rv(1:4) / sqrt(sum(rv(1:4).^2));
% rv(5) = 1;%rv(5)+0.5;
% Tnow = rv;

Fnow = mrepDist(M1, M2);

err = 1e3;
step = 5e-2;
epsilon = 1e-4; 

while (err > epsilon)
    
    gradF = transfGrad('mrepDist', M1, M2, Tnow, tMode);
    deltaT = [QuatExp(-step*gradF(1:3,:)); exp(-step*gradF(4,:))];
    Tnext = [QuatProd(deltaT(1:4,:), Tnow(1:4,:)); Tnow(5,:).*deltaT(5,:)];
    Fnext = mrepDist(applyMrepTransf(Tnext, M1), M2);
    err = abs(Fnext - Fnow);
    
    Tnow = Tnext;
    Fnow = Fnext;    
end

T1 = Tnow;

return

% =========================================================
%  Gradient method
%   
%   function [gradF] = transfGrad(F, M1, M2, T0)
%     Input: F(T; M1, M2)   - a real valued funtion 
%            T0             - point at which grad(F) is to be evaluated
%     Output: grad(F) at T0
% ---------------------------------------------------------

function [gradF] = transfGrad(F, M1, M2, T0, tMode)

step = 2e-2;
%F0 = feval(F, applyMrepTransf(T0, M1), M2);
%statFlag = 0;

T1 = T0;
T1(1:4) = QuatProd(QuatExp([step; 0; 0]), T0(1:4));
T2 = T0;
T2(1:4) = QuatProd(QuatExp([-step; 0; 0]), T0(1:4));
dx = (feval(F, applyMrepTransf(T1, M1), M2) - feval(F, applyMrepTransf(T2, M1), M2)) / step / 2;     
%dx = (feval(F, applyMrepTransf(T1, M1), M2, statFlag ) - ...
%         feval(F, applyMrepTransf(T2, M1), M2, statFlag )) / step / 2;
%dx = (feval(F, applyMrepTransf(T1, M1), M2 ) - F0) / step;
    
T1 = T0;
T1(1:4) = QuatProd(QuatExp([0; step; 0]), T0(1:4));
T2 = T0;
T2(1:4) = QuatProd(QuatExp([0; -step; 0]), T0(1:4));
dy = (feval(F, applyMrepTransf(T1, M1), M2) - feval(F, applyMrepTransf(T2, M1), M2)) / step / 2;     
%dy = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;
    
T1 = T0;
T1(1:4) = QuatProd(QuatExp([0; 0; step]), T0(1:4));
T2 = T0;
T2(1:4) = QuatProd(QuatExp([0; 0; -step;]), T0(1:4));
dz = (feval(F, applyMrepTransf(T1, M1), M2) - feval(F, applyMrepTransf(T2, M1), M2)) / step / 2;     
%dz = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;

if (tMode == 3)
    T1 = T0;
    T1(5) = T0(5) .* exp(step);
    T2 = T0;
    T2(5) = T0(5) .* exp(-step);
    dr = (feval(F, applyMrepTransf(T1, M1), M2) - feval(F, applyMrepTransf(T2, M1), M2)) / step / 2;     
%    dr = (feval(F, applyMrepTransf(T1, M1), M2, statFlag ) - ...
%             feval(F, applyMrepTransf(T2, M1), M2, statFlag)) / step / 2;
    %dr = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;
else
    dr = 0;
end
    
gradF = [dx; dy; dz; dr];

return

% % =========================================================
% %  Gradient method
% %   
% %   function [gradF] = transfGrad(F, M1, M2, T0)
% %     Input: F(T; M1, M2)   - a real valued funtion 
% %            T0             - point at which grad(F) is to be evaluated
% %     Output: grad(F) at T0
% % ---------------------------------------------------------
% 
% function [gradF] = transfGrad(F, M1, M2, T0, tMode)
% 
% step = 2e-2;
% %F0 = feval(F, applyMrepTransf(T0, M1), M2);
% 
% T1 = T0;
% T1(1:4) = QuatProd(QuatExp([step; 0; 0]), T0(1:4));
% T2 = T0;
% T2(1:4) = QuatProd(QuatExp([-step; 0; 0]), T0(1:4));
% dx = (feval(F, applyMrepTransf(T1, M1), M2 ) - ...
%         feval(F, applyMrepTransf(T2, M1), M2 )) / step / 2;
% %dx = (feval(F, applyMrepTransf(T1, M1), M2 ) - F0) / step;
%     
% T1 = T0;
% T1(1:4) = QuatProd(QuatExp([0; step; 0]), T0(1:4));
% T2 = T0;
% T2(1:4) = QuatProd(QuatExp([0; -step; 0]), T0(1:4));
% dy = (feval(F, applyMrepTransf(T1, M1), M2 ) - ...
%         feval(F, applyMrepTransf(T2, M1), M2 )) / step / 2;
% %dy = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;
%     
% T1 = T0;
% T1(1:4) = QuatProd(QuatExp([0; 0; step]), T0(1:4));
% T2 = T0;
% T2(1:4) = QuatProd(QuatExp([0; 0; -step;]), T0(1:4));
% dz = (feval(F, applyMrepTransf(T1, M1), M2 ) - ...
%         feval(F, applyMrepTransf(T2, M1), M2 )) / step / 2;
% %dz = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;
% 
% if (tMode == 3)
%     T1 = T0;
%     T1(5) = T0(5) .* exp(step);
%     T2 = T0;
%     T2(5) = T0(5) .* exp(-step);
%     dr = (feval(F, applyMrepTransf(T1, M1), M2 ) - ...
%             feval(F, applyMrepTransf(T2, M1), M2 )) / step / 2;
%     %dr = (feval(F, applyMrepTransf(T, M1), M2 ) - F0) / step;
% else
%     dr = 0;
% end
%     
% gradF = [dx; dy; dz; dr];
