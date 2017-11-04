%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SP : nSamps x 1 array of QuadSymPrimitives
% 
% atomMean(SP) computes the mean atom for all atoms in SP.
% atomMean(SP, W) computes the weighted mean atom for atoms in SP, with
%                 weights W
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function meanAtom = atomMean(SP, varargin)

nSamps = length(SP);

% convert to a vector of nSamps x 15
S = cVec(SP);

if (nargin < 2)
    W = eye(nSamps);
else
    if (length(varargin{1}) ~= nSamps)
        error('Error: weights are incorrectly assigned.');
    end
    W = diag( varargin{1} );
end

meanPos = mean(W * S(:, 1:3));

meanR = exp(mean(W * log(S(:, 4:6))));

meanU = [SphereMean(S(:,7:9)', W) ...
          SphereMean(S(:,10:12)', W) ...
          SphereMean(S(:,13:15)', W)];

meanAtom = SrepQuadPrimitive(meanPos', meanR, meanU, false);
    
return;

