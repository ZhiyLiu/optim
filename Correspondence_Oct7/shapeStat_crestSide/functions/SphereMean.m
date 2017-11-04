% [M] = SphereMean(S)
%
% Compute the mean in S^2.
% input S: (3- dim of points in S^2) x (number of points)
function [M] = SphereMean(S, W)


M = mean(S, 2);
nM = sqrt(sum(M.^2));

M = M / nM;
qM = getRotation(M);

err = 1;
epsilon = 1e-6;         % 1e-8

%h = waitbar(0, 'Computing mean...');
cntMAX = 5e4;
cnt = 0;

while (err > epsilon & cnt < cntMAX)    
    dLogS = SphereLog( QuatRotVec(qM, S) );
    dLogM = mean( dLogS*W, 2 );
    err = sum(dLogM .* dLogM);
    
    dM = SphereExp(dLogM);
    M = QuatRotVec(QuatInv(qM), dM);
    qM = getRotation(M);
    %M = QuatRotVec( QuatInv(getRotation(dM)), M );
    
    cnt = cnt + 1;
    %waitbar(max(min(1, epsilon/err), cnt/cntMAX), h, ...
    %    ['Computing mean.  Current error: ' num2str(err)]);   
end

if (cnt == cntMAX)
    disp('Warning: max iteration number achieved.');
end

%close(h); drawnow;