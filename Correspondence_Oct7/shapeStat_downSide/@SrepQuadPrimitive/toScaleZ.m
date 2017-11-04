% Use to put data to scale in z-dimension
%
% sWE: scale 

function qSrepPrim = toScaleZ(qSrepPrim, sWE)

% Scaling postion
pos = [qSrepPrim.pos(1); qSrepPrim.pos(2); sWE * qSrepPrim.pos(3)];

% Scaling normal vectors
tmpU = qSrepPrim.U;
tmpU(3,:) = sWE * tmpU(3,:);
ntmpU = sqrt(sum(tmpU.^2,1));

U = tmpU ./ repmat(ntmpU,[3,1]);

% Scaling radius
r = ntmpU .* qSrepPrim.r;

qSrepPrim = SrepQuadPrimitive(pos, r, U);

