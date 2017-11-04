% 	Input:
% 		int nObjects;
% 		int nRows;
% 		int nColumns;
% 		int nLinkAtoms;
% 		int linkAtomIds[];
% 		double oldAtoms[];
% 		double newLinkAtoms[];

% Output:
% 	    int nObjects;
% 		int nRows;
% 		int nColumns;
% 		double newAtoms[];

function [newAtoms] = predictFigure(header, oldAtoms, chgIDs, chgAtoms, dirPath)

newAtoms = [];
%tempFile = fullfile(dirPath, 'tempFig.dat');
tempFile = 'tempFig.dat';

fid = fopen(tempFile, 'wb');
if (fid == 0)
    disp('Error: Cannot open temporary file.');
    return;
end

fwrite(fid, header.nSamps, 'int');
fwrite(fid, header.nRows, 'int');
fwrite(fid, header.nCols, 'int');
fwrite(fid, length(chgIDs), 'int');
fwrite(fid, chgIDs, 'int');

[naDims, nSamps, nAtoms] = size(oldAtoms);

for k = 1:nSamps
    fwrite(fid, convert2Atom(oldAtoms(:,k,:)), 'double');
    fwrite(fid, convert2Atom(chgAtoms(:,k,:)), 'double');
end

fclose(fid);

if (ispc)
%    [status, msg] = dos(['pablo_predict.exe -pre ' tempFile], '-echo');
    searchPaths = path;
    while (1)
        [token, searchPaths] = strtok(searchPaths, ';');
        %[status, msg] = dos([token '\pablo_predict.exe -pre ' tempFile], '-echo');
        %[status, msg] = dos([token '\pablo_predict_output_sim -pretempFigNew.dat f ' tempFile], '-echo');
        [status, msg] = dos([token '\pablo_predict_output_sim -pre tempFigNew.dat f ' tempFile]);
        if(status == 0 | isempty(searchPaths))
            break;
%         elseif (status ~= 0)
%             disp('Error : pablo_predict.exe!');
%             return;
         end
    end
else
    disp('Error: Unix version of pablo_predict not available!');
    disp('Residue statistics computation not complete!');
    return;
end

% tempOutFile = fullfile(dirPath, 'tempFig_output.dat');
tempOutFile = 'tempFig_output.dat';
fid = fopen(tempOutFile, 'rb');
fseek(fid, 3*4, 'bof');

atoms = fread(fid, (naDims-1)*nSamps*nAtoms, 'double');
atoms = reshape(atoms, naDims-1, nAtoms, nSamps);
newAtoms = convert2Sym(atoms);
newAtoms = permute(newAtoms, [1 3 2]);

% otherAtomsID = setdiff(1:nAtoms, chgIDs);
% newAtoms = newAtoms(:, :, otherAtomsID);

fclose(fid);
delete(fullfile(pwd, tempFile));
delete(fullfile(pwd, tempOutFile));

return;
