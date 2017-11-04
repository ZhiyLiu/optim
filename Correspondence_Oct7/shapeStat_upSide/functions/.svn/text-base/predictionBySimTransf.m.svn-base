%  Input:
%   figIndex : prediction from figure 'figIndex' to other dependent figures
%   header : info about m-rep training model
%------------------------------------------------------------------------------
%   predictions after figure 'figIndex' is deformed
%------------------------------------------------------------------------------


function [projections, modifiedAtomIDs]  = predictionBySimTransf(header,...
    figIndex, resStruct, prevStruct, dirPath)

%check input 

numFigs = sum(header.nFigs);
projections = resStruct.projections;
modifiedAtomIDs = [];
for j = 1:numFigs
    if ( ~isempty(header.dependency.depAtomMatrix{figIndex,j}) )
        atomIDs = header.figIDs(j)+[0:header.nFigAtoms(j)-1];
        modifiedAtomIDs = [atomIDs, modifiedAtomIDs];
        depAtomIDs = header.dependency.depAtomMatrix{figIndex,j};
        
        figHeader = struct('nSamps', header.nSamps, 'nRows', header.nRows(j), ...
            'nCols', header.nCols(j));
        
         % change in atoms of the figure j
        projections(:,:,atomIDs) = predictFigure( figHeader, ...
            prevStruct.projections(:,:,atomIDs), ...
            depAtomIDs-header.figIDs(j), ...
            resStruct.projections(:,:,depAtomIDs));
        
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % propagate the prediction to other objects
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        for k = 1:numFigs
            if( ( ~isempty(header.dependency.depAtomMatrix{j,k}) ) & ( k > i ) )
                [COG, T, S, R, Q] = readSim('D:\ResidueStats-matlab-code\Functions-scripts');
                depAtomIDs = header.dependency.depAtomMatrix{j,k}; 
                depAtoms = prevStruct.projections(:, :, depAtomIDs);
                
                %----------------------------------------------
                % apply similarity transform on dependent atoms
                %----------------------------------------------
                numDepAtoms = length(depAtomIDs);
                numTargFiles = header.nSamps;
                COG = COG*ones(1, numDepAtoms*numTargFiles);
                T = T*ones(1, numDepAtoms*numTargFiles);   %numTargFiles = number of input models
                x = S*(depAtoms(1:3, :) - COG);    
                v1 = depAtoms(5:7, :);
                v2 = depAtoms(8:10, :);
                
                Rx = QuatRotVec(Q, x) + COG + T;
                Rv1 = QuatRotVec(Q, v1);
                Rv2 = QuatRotVec(Q, v2);
                Ss = S*depAtoms(4, :);     
                
                projections(1:3, :, depAtomIDs) = reshape(Rx, 3, numTargFiles, numDepAtoms);
                projections(4, :, depAtomIDs) = reshape(Ss, 1, numTargFiles, numDepAtoms);
                projections(5:7, :, depAtomIDs) = reshape(Rv1, 3, numTargFiles, numDepAtoms);
                projections(8:10, :, depAtomIDs) = reshape(Rv2, 3, numTargFiles, numDepAtoms);
                
                %----------------------------------------------
                % apply prediction to other atoms 
                %----------------------------------------------
                
                %atoms dependent on deformation of fig j
                atomIDs = header.figIDs(k)+[0:header.nFigAtoms(k)-1];                            
                modifiedAtomIDs = [atomIDs, modifiedAtomIDs];
                figHeader = struct('nSamps', header.nSamps, 'nRows', header.nRows(k), ...
                    'nCols', header.nCols(k));
                projections(:,:,atomIDs) = predictFigure( figHeader, ...
                    prevStruct.projections(:,:,atomIDs), ... %the approximation to the augmented atoms at previous scale
                    depAtomIDs-header.figIDs(k), ...
                    projections(:,:,depAtomIDs));   %new positions of the augmented atoms
%                     prevStruct.projections(:,:,depAtomIDs));                     
                
            end
        end
        
    end
end

modifiedAtomIDs = sort(modifiedAtomIDs);

return;