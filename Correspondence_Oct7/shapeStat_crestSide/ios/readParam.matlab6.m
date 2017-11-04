% =====================================================================
% Function Param = readParam(dirPath)
%
%   Read 'param.txt' file and get the setting to run 'getResidueStat.m'
%   dirPath : specifies where 'param.txt' file is stored.
% ---------------------------------------------------------------------
function Param = readParam(dirPath)

filePath = fullfile(dirPath, 'param.txt');
filename = dir(filePath);
% filename = 
% 
%      name: 'C:\MATLAB6p1\work\3106-bl-pr-re-try5-10.21.04\param.txt'
%      date: '26-Oct-2004 09:41:15'
%     bytes: 1095
%     isdir: 0

isThere = length(filename);
while (~isThere)
    beep;
    uiwait(errordlg(['No ''param.txt'' found in ' dirPath]));
    
    prompt = 'Enter path name where ''param.txt'' file is stored:';
    answer = inputdlg(prompt, 'Select directory', 1, {dirPath});
    if (isempty(answer))
        beep;
        disp('    Error reading ''param.txt'' file.  Abort.');
        return;
    else
        dirPath = answer{1};
        filePath = fullfile(dirPath, 'param.txt');
        filename = dir(filePath);
        isThere = length(filename);
    end
end


m = importdata(filePath);
strParam = strvcat( m ); % create string matrix from char array(importdata)

% read in the parameters
statsFlag = sscanf(strParam(strmatch('statsFlag', strParam), :), 'statsFlag = %d');
bReadData = sscanf(strParam(strmatch('bReadData', strParam), :), 'bReadData = %d');
PGAflag = sscanf(strParam(strmatch('PGAflag', strParam), :), 'PGAflag = %d');
bOutput = sscanf(strParam(strmatch('bOutput', strParam), :), 'bOutput = %d');
bCleanUp = sscanf(strParam(strmatch('bCleanUp', strParam), :), 'bCleanUp = %d');
prediction = sscanf(strParam(strmatch('prediction', strParam), :), 'prediction = %d');
alignMode = sscanf(strParam(strmatch('alignMode', strParam), :), 'alignMode = %d');
normalizeEndAtoms = sscanf(strParam(strmatch('normalizeEndAtoms', strParam), :), 'normalizeEndAtoms = %d');

strFigOrder = strParam(strmatch('dependency.figureOrder', strParam), :);
[tok1, strFigOrder] = strtok(strFigOrder, '[');
strFigOrder = strFigOrder(1, 2:length(strFigOrder));
figureOrder = [];
while (1)
    [tok1, strFigOrder] = strtok(strFigOrder, ' ]');
    if (isempty(tok1)), break; end  
    a1 = str2num(tok1);
    if (isempty(a1)), continue; end    
    figureOrder = [figureOrder, a1];
end
numFigs = length(figureOrder);

strFigNames = strParam(strmatch('figNames', strParam), :);
[tok4, strFigNames] = strtok(strFigNames, '[');
strFigNames = strFigNames(1, 2:length(strFigNames));
figNames = cell(1, numFigs);
h = 1;
while (1)
    [tok4, strFigNames] = strtok(strFigNames, ',] ');
    if (isempty(tok4) | strfind(tok4, ';')), break; end
    figNames{h} = tok4;
    h = h+1;
end
figNames = {figNames'};

if (h-1 ~= numFigs & numFigs > 0)  %if numFigs == 0, only global statistics
    errordlg('Parameters are not consistent. Check ''param.txt'' !');
    return;
end

rem = strParam;
[nr, nc] = size(rem);
rem = reshape(rem', 1, nr*nc);
[tok, rem] = strtok(rem, '{');
[tok, rem] = strtok(rem, ';');
strdepAtomTree = deblank(tok);
i = 1;
depAtomTree = cell(1, numFigs);
while (1)
    
    [tok, strdepAtomTree] = strtok(strdepAtomTree, ',');
    if (isempty(tok)), break; end 
    tok = deblank(tok);
    
    remFig = tok;
    j = 1;
    while (1)
        depAtom = [];
        [tok, remFig] = strtok(remFig, '|');
        if (isempty(tok)), break; end
        
        %remove {
        rem = tok;
        while (1) 
            [tok, rem] = strtok(rem, '{');
            if (isempty(rem)), break; end
        end
        
        subTok = tok;
        while (1) 
            [tok, subTok] = strtok(subTok, ' }');
            if (isempty(subTok)), break; end
            depAtom = [depAtom str2num(tok)];
        end
        %depAtom
        %i, j
        %depAtomTree{i}{j} = mat2cell(depAtom);  %single input to mat2cell will be
        %obsolete in future release of matlab
        if (length(depAtom))
            depAtomTree{i}{j} = depAtom;
        end
        j = j+1;
    end
    i = i+1;
end

depAtomTree = {depAtomTree};

% if (i-1 ~= numFigs & numFigs > 0) % if numFigs == 0, only global statistics
%     errordlg('Parameters are not consistent. Check ''param.txt'' !');
%     return;
% end
    
strnComp = strParam(strmatch('nComp', strParam), :);
[tok2, strnComp] = strtok(strnComp, '[');
strnComp = strnComp(1, 2:length(strnComp));
nComp = [];
while (1)
    [tok2, strnComp] = strtok(strnComp, ' ]');
    if (isempty(tok2)), break; end  
    a2 = str2num(tok2);
    if (isempty(a2)), continue; end    
    nComp = [nComp, a2];
end

strnProj = strParam(strmatch('nProj', strParam), :);
% JEONG : this is weird
%strnProj = strtok(strnProj, ';')
strnProj = strnProj(1, :);

[tok3, strnProj] = strtok(strnProj, '[');
strnProj = strnProj(1, 2:length(strnProj));
nProj = [];
while (1)
    [tok3, strnProj] = strtok(strnProj, ' ]');
    if (isempty(tok3)), break; end  
    a3 = str2num(tok3);
    if (isempty(a3)), continue; end    
    nProj = [nProj, a3];
end

if ( ~length(statsFlag) | ~length(bReadData) | ~length(PGAflag) | ~length(bOutput) | ... 
     ~length(bCleanUp) | ~length(nComp) | ~length(nProj)  | ~length(figNames))
    errordlg('All parameters are not specified. Check ''param.txt'' !');
    return;
end

%Param = struct('statsFlag', statsFlag, 'bReadData', bReadData, 'PGAflag', PGAflag, 'bOutput', bOutput,...
%    'bCleanUp', bCleanUp, 'figNames', figNames, 'figureOrder', figureOrder, 'depAtomTree', depAtomTree,  'nComp', nComp, 'nProj', nProj);

Param = struct('statsFlag', statsFlag, 'bReadData', bReadData, 'PGAflag', PGAflag, 'bOutput', bOutput, ...
    'bCleanUp', bCleanUp, 'figNames', figNames, 'figureOrder', figureOrder, 'depAtomTree', depAtomTree, ... 
    'nComp', nComp, 'nProj', nProj, 'prediction', prediction, 'alignMode', alignMode,
	'normalizeEndAtoms', normalizeEndAtoms );

% ===========  Set these parameters to appropriate values  ====================

% statsFlag = 1;          % 0-Lie group, 1-symmetric space.  DON'T CHANGE THIS!
% 
% bReadData = 0;          % 0-skip, 1-read
% PGAflag = 1;            % -1-skip, 0-all (no atom!!), 1-ensemble  
%                         % Change dependency and nComp & nProj accordingly
%                         
% bOutput = 2;            % 0-skip, 1-binary output, 2-ascii output
% bCleanUp = 1;           % 0-no clean up, 1-clean up


	% ----------------------------------------------------------------------------
	% Dependency graph  
	%
	%   dependency is a structure with the following fields:
	%     - figureOrder:    each entry is a figure number.  The statistics are 
	%                       computed in the order given here.
	%     - depAtomTree:    a cell array with nested cell entries.  The first
	%                       level is an array with the same length as figureOrder.
	%                       Each entry is itself a cell array.
	%
	% 
	% if PGAflag = 0, i.e., in case getting statistics for ensemble stage and 
	% each objects with augmented atoms,
	%   for two objects(figure number 0, 1) and atoms 1, 2, 5 of figure 1 are augmented 
	%   to object of figure 0, 
	%   dependency.figureOrder = [ 0 1 ];
	%   dependency.depAtomTree = {{{1,1,2,5}}, ... //{figure 1, atom1, atom2, atom 5}
	%                              {}, ... 
	%                             };
	%   for four objects(figure number 0, 1, 2, 3) with augmented atoms
	%   dependency.figureOrder = [0 1 2 3 ];
	%   dependency.depAtomTree = {{{1, 2,3,4}},...
	%                             {{2, 3,5,7}},...
	%                             {{3, 2,3,6}},...
	%                             {}};
	% ----------------------------------------------------------------------------

% (global, fig1, fig2)
% nComp = [8];            % number of PG's to be computed for each PGA 
% nProj = [3];            % number of PG's used for computing residues.  For 
                         % example, if nProj=5, the first 5 PG's are taken out
                        
% ====================  End of setting parameters  ===========================
