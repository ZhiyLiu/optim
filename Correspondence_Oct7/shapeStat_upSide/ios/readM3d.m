function [header, s] = readM3d(m3dDir, filetype, varargin)
% Read a set of .m3d files from a specified directory.  
%
%    [header, s] = readM3d(m3dDir, filetype,  varargin)
%
% Returns:
%   header: Header information
%   s: Information about primitives and other things.
%
% Parameters: 
%   m3dDir: The directory containing the files
%   filetype: What files to read in that directory, described using 
%        wildcards.  Usually '*.m3d'.  Can just be a filename.
%   varargin: An optional 3rd argument is a struct of options.  Currently
%        they are:
%          * 'list': a list of indices specifying which of the filenames 
%             found in the directory to actually read.
%          * 'reread': If true, re-read the data from the m3d files rather
%             than from the M3D.dat file (default is False)
%          * 'makeM3dDataFile: If true, write out the data that is loaded
%            to M3D.dat.  (Default is True; make false to leave existing
%            M3D.dat untouched.)
% 
% Example: 
%  To read model files of the form test/models/run1*.m3d, but only the 1st,
%  2nd, and 7th files of that form in the directory, do this:
%   opt.list = [1 2 7];
%   [h, s] = readM3d('test/models', 'run1*.m3d', opt.list);

%% Some global constants (enums for tube and quad figures)
% The values correspond to those in the cell array.
%
global QUAD_FIGURE;
global TUBE_FIGURE;
global FIG_TYPENAMES;
QUAD_FIGURE = 1;
TUBE_FIGURE = 2;
FIG_TYPENAMES   = { 'QuadFigure', 'TubeFigure' };

reread = true; %% always read from m3d file.
makeM3dDataFile = true;
if (nargin > 2)
    options = varargin{1};
    if isfield(options, 'list')
        fileList = options.list;
    end
    if isfield(options, 'reread')
        reread = options.reread;
    end
    if isfield(options, 'makeM3dDataFile')
        makeM3dDataFile = options.makeM3dDataFile;
    end
end


% file to store point list of m3d tiles
m3dDatafile = fullfile(m3dDir, 'M3D.mat');

%% If m3d files are already read, then load the previously read data.
if ( ~reread && exist(m3dDatafile, 'file') )
    % Check the number of files in the input data directory (m3dFiles) and
    % the number of files already stored in 's'. If they are not the same,
    % reread the m3d files in the input data directory.
    m3dFiles = fullfile(m3dDir, filetype);
    fStructs = dir(m3dFiles);
    nFiles = size(fStructs, 1);
    
    disp('   load previously read data of *.m3d files ... ');
    disp(' ');
    load(m3dDatafile, 'header', 's');
    if (nFiles ~= size(s.atoms(:,1), 1))
        disp('  Previously read data is not consistent with the current data ... reloading files.');
        disp('');
    else
        return;
    end
end

%--------------------------------------------------------------------
% If not, read.
%--------------------------------------------------------------------

m3dFiles = fullfile(m3dDir, filetype);
pathstr = fileparts(m3dFiles);
fStructs = dir(m3dFiles);

if (exist('fileList', 'var'))
    fStructs = fStructs(fileList);
end

nFiles = size(fStructs, 1);

s = [];

while (nFiles == 0)
    beep;
    uiwait(errordlg(['No m3d files found in ' pathstr]));
    
    prompt = 'Enter path name where m3d files are stored:';
    answer = inputdlg(prompt, 'Select directory', 1, {m3dFiles});
    charAnsw = char(answer);
    if (isempty(answer))
        beep;
        disp('    Error reading m3d files.  Abort.');
        return;
    elseif (~strcmp(charAnsw((length(charAnsw)-4):length(charAnsw)), '*.m3d'))
        beep;
        disp('    Please add ''*.m3d'' at the end of the path to the data directory.');
        disp(' ');
    else
        m3dFiles = answer{1};
        pathstr = fileparts(m3dFiles);
        fStructs = dir(m3dFiles);
        nFiles = size(fStructs, 1);
    end
end

disp(['    ' num2str(nFiles) ' m3d file(s) found.']);

% default: arg is (m3dFiles)
fileSpec = '';

m3dCellArray = [];
simTrans = [];
elVec = [];
worldExts = cell(1, nFiles);
maxDim = []; % store the maximum extent among x, y, z and also used as flag whether to
             % take the world extent into account in computing statistics.

fNames = cell(1, nFiles);
h = waitbar(0, ['Reading ' fileSpec '...']);

%% Read Header information from the first file 

% Note : Does NOT parse multi-figure object m3d file correctly.
fNames{1} = fullfile( pathstr, fStructs(1).name );
landmarks = readLandmarks(fNames{1});
figNames = readFigNames(fNames{1});  

% Use flat_d.exe to keep all the extra parameters in the mean file 
flatCArrayfile1 = flat(fNames{1});

% create char array % fix for version 7.0
% Ja-Yeon 3/30/2005 delimeter '\n', '\t' in importdata()....
sTest = strvcat( importdata(fNames{1}, '\n') );    
stemp = removeBlanks(sTest);
%NOTE: Lines below only work until Matlab ver.6.5.  
%importdata(), sprintf() have changed in ver. 7.x
%stemp = strvcat( importdata(fNames{i}) );  % create char array 

% fix for version 7.0
nChilds = readLines('childCount', stemp)';
nObjs = readLines('count', stemp)';
nCols = readLines('numColumns', stemp)';
nRows = readLines('numRows', stemp)';

figStringTypes  = readStringLines('type', stemp, 'Figure');
figTypes    = zeros([1, length(figStringTypes)]);
figTypes( strcmp( figStringTypes, FIG_TYPENAMES{QUAD_FIGURE} ) )    = QUAD_FIGURE;
figTypes( strcmp( figStringTypes, FIG_TYPENAMES{TUBE_FIGURE} ) )    = TUBE_FIGURE;

[parentFigIds rootFigIds] = readTree(flatCArrayfile1, nObjs, nChilds);

nFigs = zeros(1, nObjs);
for i = 1:nObjs
    nFigs(i) = sum(nChilds(find(rootFigIds == i-1)))+1;
end


nFigAtoms = nCols .* nRows;     % number of atoms in each figure
figIDs = cumsum( [0 nFigAtoms(1:length(nFigAtoms) - 1)] ) + 1;  % first atom in each figure

% This is the host figure ID for each object
objIDs = cumsum([0 nFigs(1:nObjs-1)]) + 1;

colors = [readLines('blue', stemp), ...
          readLines('green', stemp), ...
          readLines('red', stemp) ];
%------- End of Reading Header information---------------------------------------

%------- Reading primitives -----------------------------------------------------
noWorldExts = false;
for i = 1:nFiles
    
    fNames{i} = fullfile( pathstr, fStructs(i).name );
    stemp = flat(fNames{i});
    
    %% Read atom information
    % parameters of atom : ((x,y,z) r, (qw, qx, qy, qz), theta, elongation) 10
    % parameters!
    
    prims = parsePrimLines(stemp, nFigs, nCols, nRows);
    simTrans = [simTrans parseSimTrans(stemp)];
    
    % JJ: To deal with the need to access atoms of one case in the training sample
    m3dCellArray = [ m3dCellArray; prims ];  		
       
    if (~noWorldExts)
        worldExts{i} = parseExt(stemp);
        if (isempty(worldExts{i})) % if one of training models does not have world extent, 
                                   % do not take world extent into account in computing statistics.
            noWorldExts = true;
        elseif (~noWorldExts)
            maxDim = [maxDim; max(worldExts{i}{6})];
        end
    else
        worldExts{i} = {};
        maxDim = [];
    end
    
    waitbar(i/nFiles, h, ['Reading ' fileSpec ': ' ...
            num2str(i) ' out of ' num2str(nFiles) ' done...']);
end

close(h);
freqDim = [];
if (~noWorldExts)
    freqDim = findMostFreqDim(maxDim); 
    disp(' ');
    disp(['    Most frequent maximum extent in models is ' num2str(freqDim(1)) '.']);
    disp(['    Selected model file for world extents: ' fNames{freqDim(2)} '.']);
end


% ------------------------------------------------------------------------
% Format of m3dArray: cell array of column vectors of each atom
%      A1 A2 ... An
%      B1 B2 ... Bn
%          ...
%      M1 M2 ... Mn
% ------------------------------------------------------------------------
[sampleSize, nAtoms] = size(m3dCellArray);
% for i = 1: nAtoms
%     aPrimColumn = [];
%     for j = 1: sampleSize
%         aPrimColumn = [ aPrimColumn ; m3dCellArray{j,i} ];
%     end
%     m3dCellArrayColumn(i) = {aPrimColumn};
% end


%% Create struct
header = struct('nObjs', nObjs, 'nFigs', nFigs, 'nAtoms', nAtoms, ...
	'figTypes', figTypes, 'nChilds', nChilds, ...
    'nRows', nRows, 'nCols', nCols, 'nSamps', nFiles, 'colors', colors, ...
    'nFigAtoms', nFigAtoms, 'objIDs', objIDs, 'figIDs', figIDs, ...
    'landmarks', landmarks, 'figNames', figNames, 'flatCArrayfile1', flatCArrayfile1);
s = struct('atoms', {m3dCellArray}, 'simTrans', simTrans, 'rawFileNames', {fNames}, ...
    'worldExtents', {worldExts}, 'maxDim', maxDim, 'freqDim', freqDim);

%--------------------------------------------------------------------
% Save read data
%--------------------------------------------------------------------

if makeM3dDataFile
    save(m3dDatafile, 'header', 's');
    disp(' ');
    disp([ '   M-reps in ' m3dFiles ' are stored at ' m3dDatafile '.' ]);
end

return;
