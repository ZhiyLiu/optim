function [ predStruct, projections, remAtomIDs] = predictionByShapeSpace(workingDir, header, figIndex,...
                                                  atomIDs,...
                                                  prevStruct, resStruct, ...
                                                  nComp, nProj, statsFlag, meanRadii)

numFigs = sum(header.nFigs);

if ( figIndex+1 > numFigs)
    disp('    Error : last Object. Nothing to predict!');   
    return;
else
    j = header.dependency.figureOrder(figIndex+1);
    remAtomIDs = [];
    nRemAtoms = 0;
    for k = figIndex+1:numFigs
        idx = header.dependency.figureOrder(k);
        remAtomIDs = [remAtomIDs, header.figIDs(k) + [0:header.nFigAtoms(k)-1]];
        nRemAtoms = nRemAtoms + header.nFigAtoms(k);
    end
    remAtomIDs = sort(remAtomIDs);    
   
    if(length(remAtomIDs) ~= nRemAtoms)
        disp('    Error : mismatch in atoms IDs and number of atoms in the remainder!');                
    end
    
    fileName = fullfile(workingDir,  ['eModes_Pred_from_Obj' num2str(figIndex) '.jpg']);
    tempStruct = getPGA(prevStruct.residues(:,:,remAtomIDs), nComp, nProj, ...
        statsFlag, fileName, meanRadii(remAtomIDs));
    
    %-------------------------------------------------------------
    % Store statistics for prediction in 'predStruct(ii)'
    %-------------------------------------------------------------
    predStruct = prevStruct;
    predStruct.stats.Mean(:,remAtomIDs) = tempStruct.stats.Mean;
    predStruct.stats.PCs(:,remAtomIDs,:) = tempStruct.stats.PCs;
    predStruct.stats.EVs = tempStruct.stats.EVs;
    predStruct.logProjs = tempStruct.logProjs;
    predStruct.numPGs = tempStruct.numPGs;
    predStruct.projDims = tempStruct.projDims;
    
    %-------------------------------------------------------------
    % compute the prediction
    %-------------------------------------------------------------
    depAtomIDs = [header.dependency.depAtomMatrix{figIndex,:}];
    projections = resStruct.projections;
    projections(:,:,remAtomIDs) = predictCoeff( atomIDs, depAtomIDs, remAtomIDs, ...
        resStruct, ... 
        predStruct, ...
        prevStruct.projections);   
    
end