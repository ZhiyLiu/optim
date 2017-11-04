% =====================================================
%  Gradient descent method
%   
%   function {Qopt] = rotGradDescent(F, Q0, aMatrix, meanR)
%     Input: F       - real valued funtion to be optimized
%            Q0      - initial starting point
%            aMatrix - atom matrix
%            meanR   - mean radii
%     Output: the optimizer
% -----------------------------------------------------

function [Q] = rotGradDescent(F, Q0, aMatrix, meanR)

Qnow = Q0;
fnow = feval(F, Qnow, aMatrix, meanR);

err = 1e3;
step = 1e-2;
epsilon = 1e-4;

%h = msgbox('Aligning models.  Please wait...', 'Aligning models');
h = waitbar(0, 'Aligning models.  Please wait...');

while (err > epsilon)
    gradF = rotGrad(F, Qnow, aMatrix, meanR);
    deltaQ = [QuatExp(-step*gradF(1:3,:)); max(1, exp(-step*gradF(4,:)))];
    Qnext = [QuatProd(Qnow(1:4,:), deltaQ(1:4,:)); Qnow(5,:).*deltaQ(5,:)];
    fnext = feval(F, Qnext, aMatrix, meanR);
    err = abs(fnext - fnow);
    waitbar(min(1, epsilon/err), h, ['Current error: ' num2str(err)]);
    
    Qnow = Qnext;
    fnow = fnext;
end

close(h);

Q = Qnow;

% =========================================================
%  Gradient method
%   
%   function {gradQ] = rotGrad(F, Q0, aMatrix, meanR)
%     Input: F       - a real valued funtion 
%            Q0      - point at which grad(F) is to be evaluated
%            aMatrix - atom matrix
%            meanR   - mean radii
%     Output: grad(F) at Q0
% ---------------------------------------------------------

function [gradQ] = rotGrad(F, Q0, aMatrix, meanR)

step = 5e-2;
F0 = feval(F, Q0, aMatrix, meanR);
gradQ = [];

for i = 1:9
    Q = Q0;
    Q(1:4, i) = QuatProd(Q0(1:4, i), QuatExp([step 0 0]'));
    dx = (feval(F, Q, aMatrix, meanR ) - F0) / step;
    
    Q = Q0;
    Q(1:4, i) = QuatProd(Q0(1:4, i), QuatExp([0 step 0]'));
    dy = (feval(F, Q, aMatrix, meanR ) - F0) / step;
    
    Q = Q0;
    Q(1:4, i) = QuatProd(Q0(1:4, i), QuatExp([0 0 step]'));
    dz = (feval(F, Q, aMatrix, meanR ) - F0) / step;

    Q = Q0;
    Q(5, i) = Q0(5, i) .* exp(step);
    dr = (feval(F, Q, aMatrix, meanR) - F0) / step;
    
    gradQ = [gradQ [dx; dy; dz; dr]];
end
