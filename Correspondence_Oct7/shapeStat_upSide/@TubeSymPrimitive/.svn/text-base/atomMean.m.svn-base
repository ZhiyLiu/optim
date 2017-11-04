%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SP : nSamps x 1 array of TubeSymPrimitives
% 
% atomMean(SP) computes the mean atom for all atoms in SP.
% atomMean(SP, W) computes the weighted mean atom for atoms in SP, with
%                 weights W
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function meanAtom = atomMean(SP, varargin)

nSamps = length(SP);
baseAtom	= SP(1).baseAtom;

% convert to a vector of nSamps x 9
S = cVec(SP);

if (nargin < 2)
    W = eye(nSamps);
else
    if (length(varargin{1}) ~= nSamps)
        error('Error: weights is incorrectly assigned.');
        return;
    end
    W = diag( varargin{1} );
end

meanDr	= mean([SP.dr] .* repmat(diag(W)', [size(SP(1).dr,1), 1]), 2 );

MA = [	mean(S(:,1).*diag(W)); ...
		mean(S(:,2).*diag(W)); ...
		mean(S(:,3).*diag(W)); ...
        exp( mean(log(S(:,4)).*diag(W)) ); ...
		SphereMean(S(:,5:7)', W); ...
		atan(mean(tan(S(:,8) - pi/2.0).*diag(W)) ) + pi/2.0; ...
        exp( mean(log(S(:,9)).*diag(W)) ) ];
 
%(pos: 3x1 vector, r, elongation, U0: 3x1(or 2x1) vector,
%                          InTangentSpace: boolean
meanAtom = TubeSymPrimitive( MA(1:3), MA(4), MA(9), MA(5:7), MA(8), baseAtom, meanDr, false);
 
return;

