function writeM3d(fileName, header, modelType, aMatrix, varargin)
% Create m3d model file from list of atoms
% 
%   writeM3d(fileName, header, modelType, aMatrix, varargin)
%
%   fileName:  Name of file to save to
%   header:    Header as returned by readM3d
%   modelType: Must be -1, 0, or 1.  Don't know what it means.  -1 is the
%       default.
%   aMatrix: Cell array of tube and quad primitives [1 x #(atoms)].  If the
%       model was loaded by [h, s] = readm3d(...), then aMatrix should be
%       s.atoms.
%   A final optional parameter can be either an array of indices of
%   primitives to write, or a cell array specifying the world dimensions
%   (see worldExt below for the format of this).

% mStruct.atoms
nFigs = header.nFigs;
nCols = header.nCols;
nRows = header.nRows;
figIDs = header.figIDs;
nFigAtoms = header.nFigAtoms;

if (modelType == -1) % when pga stat is not necessary
    rows = strmatch('PGAStats', header.flatCArrayfile1);
    if (isempty(rows))
        flatStrArray = cellstr(header.flatCArrayfile1);
    else
        flatStrArray = cellstr(header.flatCArrayfile1(1:min(rows)-1, :));
    end
else
    flatStrArray = cellstr(header.flatCArrayfile1);
end

[tmp, modelName, ext] = fileparts(fileName);
flatStrArray = regexprep(flatStrArray, 'model\.name .+$', ...
    ['model.name ' modelName]);
flatStrArray = regexprep(flatStrArray, 'figure\[([0-9]+)\]\.name .+$', ...
    ['figure[$1].name ' modelName]);

worldExt = {};
if (nargin == 5)
    if (isnumeric(varargin{1}))
        aMatrix = squeeze(aMatrix(varargin{1}, :));
    else
        worldExt = varargin{1};
    end
end

if (~ismember(modelType, -1:1))
    modelType = -1; %default
end

flatStrArray = [flatStrArray ; cellstr(['model.modelType ' num2str(modelType)])];

for j = 1:sum(nFigs)
    baseStr = ['model.figure[' num2str(j-1) '].'];
    aMatrixPerFig = aMatrix([figIDs(j):figIDs(j)+nFigAtoms(j)-1]);
	flatStrArray = [flatStrArray; write( reshape([aMatrixPerFig{:}], [nCols(j), nRows(j)])', baseStr ) ];
end


if (~isempty(worldExt))
	worldBase = 'model.world.';
	flatStrArray = [flatStrArray ; cellstr([worldBase  'imageModTime ' num2str(worldExt{1})])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'imagePath ' num2str(worldExt{2})])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'bound.x ' num2str(worldExt{3}(1))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'bound.y ' num2str(worldExt{3}(2))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'bound.z ' num2str(worldExt{3}(3))])];

	flatStrArray = [flatStrArray ; cellstr([worldBase  'origin.x ' num2str(worldExt{4}(1))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'origin.y ' num2str(worldExt{4}(2))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'origin.z ' num2str(worldExt{4}(3))])];

	flatStrArray = [flatStrArray ; cellstr([worldBase  'spacing.x ' num2str(worldExt{5}(1))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'spacing.y ' num2str(worldExt{5}(2))])];
	flatStrArray = [flatStrArray ; cellstr([worldBase  'spacing.z ' num2str(worldExt{5}(3))])];
end


% ************************************
% Generate the m3d model file.       *
% ************************************

nRows = length(flatStrArray);

flatfilename = [fileName, '.txt'];

fid = fopen(flatfilename, 'wt');

if (fid == -1)
    disp(['    Error: cannot open file ' fileName '.']);
else
    for i=1:nRows
        fprintf(fid, '%s\n', flatStrArray{i});
    end
    fclose(fid);
end

%fclose(fid);

if (ispc)
    command = ['flat_d.exe -u ' fileName ' ' flatfilename];
    searchPaths = path;
    while (1)
        [token, searchPaths] = strtok(searchPaths, ';');
        [status, msg] = dos([token '\' command]);
        if(status == 0 || isempty(searchPaths))
            fprintf(1, '    Output %s\n', fileName);
            break;
        end
    end
elseif (isunix)
    command = ['linflat -u ' fileName ' ' flatfilename];
    searchPaths = path;
    while (1)
        [token, searchPaths] = strtok(searchPaths, ':');
        [status, msg] = unix([token '/' command]);
        if (status == 0 || isempty(searchPaths))
            fprintf(1, '    Output %s\n', fileName);
            break;
        end
    end
else
    disp('Error: This is not a PC or Unix compatible machine!');
    disp('Returning from writeM3d.m. Residue statistics computation incomplete!');
    return;
end

% delete(flatfilename);


