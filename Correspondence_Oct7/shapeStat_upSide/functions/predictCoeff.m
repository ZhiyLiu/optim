% =======================================================================================================
% Function
%   projections = predictCoeff(atomIDs, depAtomIDs, remAtomIDs, resStruct, predStruct, prevProjs)
%  
%   Input:
%       - atomIDs:     atom IDs of the base object from which the  prediction is made
%       - depAtomIDs:  augmented atom IDs 
%       - remAtomIDs:  atom IDs in the remainder
%       - resStruct:   stats of base atoms
%       - predStruct:  stats of remainder atoms
%       - prevProjs:   best approximation to the remainder(all objects o/t base object) before prediction 
%   Output:
%       - projections: best approximation to the remainder(all objects o/t base object) after prediction
% --------------------------------------------------------------------------------------------------------

function projections = predictCoeff(atomIDs, depAtomIDs, remAtomIDs, resStruct, predStruct, prevProjs)

[nlDims, bnSamps, bnAtoms] = size(resStruct.logProjs);
[nlDims, tnSamps, tnAtoms] = size(predStruct.logProjs);

if ( bnSamps ~= tnSamps )
    disp(' error in predictFigureNew');
    disp('\n');
    return;
end
naDims = nlDims + 2;
nSamps = bnSamps;
nDepAtoms = length(depAtomIDs);
%nbAtoms = length(atomIDs);

bnComps = resStruct.projDims;
tnComps = predStruct.projDims;
if (bnComps < tnComps)
    nComps = bnComps;
else
    nComps = tnComps;
end


bMean = resStruct.stats.Mean(:, depAtomIDs);
bAugAtoms = resStruct.projections(:,:,depAtomIDs);

%bLogMean = repmat(atomLog(bMean), [1, 1, nSamps]);
tMean = predStruct.stats.Mean(:, depAtomIDs);
tMean = reshape(repmat(tMean, nSamps, 1), [naDims, nSamps, nDepAtoms]);  
tAugAtomRes = atomDiff(tMean, bAugAtoms); %naDims, nSamps, nAtoms

tLogAugAtomRes = atomLog(tAugAtomRes);
tLogAugAtomRes = reshape(tLogAugAtomRes, [nlDims, nSamps, nDepAtoms]);

tPC = predStruct.stats.PCs(:, depAtomIDs, [1:nComps]);
tEV = predStruct.stats.EVs(1:nComps);
tPC = reshape(tPC, nlDims*nDepAtoms, nComps)./ repmat(sqrt(tEV), nlDims*nDepAtoms, 1);

tLogAugRes = reshape(permute(tLogAugAtomRes, [1 3 2]), nlDims*nDepAtoms, nSamps);
%tPComp = tPC'*tLogAugRes;
newPComp = tPC'*tLogAugRes;  %nComp X nSamps

% *******************************************************************************   
% Use predicted projection('newProj') 
% *******************************************************************************   
% PC : nlDims*tnAtoms X nComps, newProj : nComps X nSamps
PCs = reshape(predStruct.stats.PCs(:, remAtomIDs, [1:nComps]), nlDims*tnAtoms, nComps); 
aLogPGProj = reshape( PCs * (newPComp .* repmat(sqrt(tEV'), 1, nSamps)), nlDims, tnAtoms*nSamps );
aPGProj = atomExp(aLogPGProj);
naDims = nlDims + 2;
aPGProj = permute( reshape(aPGProj, naDims, tnAtoms, nSamps), [1 3 2] );

projMean = atomMult(prevProjs(:,:,remAtomIDs), ...
            reshape(repmat(predStruct.stats.Mean(:, remAtomIDs), nSamps,1), ... *nlDim*nSamps X nAtoms 
                    size(prevProjs(:,:,remAtomIDs))));
projections = atomMult(projMean, aPGProj);    
        
% check the units!!!
% how to check the correctness of computation done here?
% adjustUnit??




                               
                              