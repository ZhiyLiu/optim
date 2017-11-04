% Output:
% 	    int nObjects;
% 		int nRows;
% 		int nColumns;
% 		double newAtoms[];

function [newAtoms] = updateFigure(oldAtoms, chgIDs, chgAtoms)

[naDims, nSamps, nAtoms] = size(oldAtoms);

fid = fopen('figDataAtom_output.dat', 'rb');
fseek(fid, 3*4, 'bof');

atoms = fread(fid, (naDims-1)*nSamps*nAtoms, 'double');
atoms = reshape(atoms, naDims-1, nAtoms, nSamps);
newAtoms = convert2Sym(atoms);
newAtoms = permute(newAtoms, [1 3 2]);

% otherAtomsID = setdiff(1:nAtoms, chgIDs);
% newAtoms = newAtoms(:, :, otherAtomsID);

fclose(fid);

return;

% ----------------------------------------------------------------------------------------
% oldCOG = mean(oldAtoms(7:9, :, :), 3);
% trans = mean( chgAtoms(7:9, :, :) - oldAtoms(7:9, :, chgIDs), 3 );
% newCOG = oldCOG + trans;
% 
% oldAtoms(7:9, :, :) = oldAtoms(7:9, :, :) - repmat(oldCOG, [1 1 nAtoms]);
% chgAtoms(7:9, :, :) = chgAtoms(7:9, :, :) - repmat(newCOG, [1 1 length(chgIDs)]);
% 
% T = [];
% for i = 1:nSamps
%     T = [T alignMreps(oldAtoms(:,i,chgIDs), chgAtoms(:,i,:), 3)];
% end
% 
% otherAtomsID = setdiff(1:nAtoms, chgIDs);
% newAtoms = oldAtoms(:, :, otherAtomsID);
% for i = 1:nSamps
%     newAtoms(:, i, :) = applyMrepTransf(T(:, i), squeeze(oldAtoms(:, i, otherAtomsID)));
% end
% newAtoms(7:9, :, :) = newAtoms(7:9, :, :) + repmat(newCOG, [1 1 length(otherAtomsID)]);
