%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SP : nSamps x 1 array of QuadSymPrimitives
% 
% atomMean(SP) computes the mean atom for all atoms in SP.
% atomMean(SP, W) computes the weighted mean atom for atoms in SP, with
%                 weights W
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function meanAtom = atomMean(SP, varargin)

nSamps = length(SP);

% convert to a vector of nSamps x 11
S = cVec(SP);

if (nargin < 2)
    W = eye(nSamps);
else
    if (length(varargin{1}) ~= nSamps)
        error('Error: weights is incorrectly assigned.');
    end
    W = diag( varargin{1} );
end

MA = [mean(S(:,1).*diag(W));  mean(S(:,2).*diag(W)); mean(S(:,3).*diag(W)); ...
        exp( mean(log(S(:,4)).*diag(W)) ); SphereMean(S(:,5:7)', W); ...
        SphereMean(S(:,8:10)', W); exp( mean(log(S(:,11)).*diag(W)) ) ];
    
%(pos: 3x1 vector, r, elongation, Up1: 3x1(or 2x1) vector,
%                          Um1: 3x1(or 2x1) vector, inTangentSpace: boolean
meanAtom = QuadSymPrimitive( MA(1:3), MA(4), MA(11), MA(5:7), MA(8:10), false);
    
return;

