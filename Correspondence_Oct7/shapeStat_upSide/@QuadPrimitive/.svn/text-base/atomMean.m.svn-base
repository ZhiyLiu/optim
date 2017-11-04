%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% P : nSamps x 1 array of QuadPrimitives
% 
% atomMean(P) computes the mean atom for all atoms in SP.
% atomMean(P, W) computes the weighted mean atom for atoms in SP, with
%                 weights W
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


function meanAtom = atomMean(P, varargin)

nSamps = length(P);

% convert to vector of nSamps x 10
S = cVec(P);  

% Get weights.  Default is equal weights:
if (nargin < 2)
    W = eye(nSamps);
else
    if (length(varargin{1}) ~= nSamps)
        error('Error: weights is incorrectly assigned.');
        return;
    end
    W = diag( varargin{1} );
end

% Compute mean rotation
M = S(1, 5:8)';       %initial guess
err = 1;
epsilon = 1e-4;     %1e-8

while (err > epsilon)
    dLogS = QuatLog( QuatProd(QuatInv(M), S(:, 5:8)') );
    dLogM = mean( dLogS*W, 2);
    dM = QuatExp(dLogM);
    err = sum(sum(dLogM .* dLogM));
    
    M = QuatProd(M, dM);
end

meanAtom = QuadPrimitive([mean(S(:,1).*diag(W)); mean(S(:,2).*diag(W)); mean(S(:,3).*diag(W))], ...
    exp( mean(log(S(:,4)).*diag(W)) ), exp( mean(log(S(:,10)).*diag(W) ) ), ...
    M, mean(S(:,9).*diag(W)));
           
