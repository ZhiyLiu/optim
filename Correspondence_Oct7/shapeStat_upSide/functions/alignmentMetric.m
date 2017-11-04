% ======================================================================
% FUNCTION [d] = alignmentMetric(mrepMatrix)
% 
%   Compute the alignment metric among mreps in mrepMatrix.
%
% INPUT: mrepMatrix - cell array of combined primitives (tube + quad)
% ----------------------------------------------------------------------

function [d] = alignmentMetric(mrepMatrix)

[nObjs, nAtoms] = size(mrepMatrix);
% meanR = [];
% for i = 1:nAtoms
%     meanR = [meanR geomean( mrepMatrix(4,:,i) )];
% end
d = 0;

%jeong 11/5/04
%statFlag = 1;

for i = 1:nObjs-1
    for j = i+1:nObjs
% jeong 11/5/04 : meanR vs statFlag  
%                 meanR to commensurate units? 
         d = d + mrepDist(mrepMatrix(i,:), mrepMatrix(j,:));
    end
end


