% ===========================================================================================
% Write PGA into a file. Note: only for symmetric representation.
%
%  - writePGA(fileName, header, PGAflag, type, globalStruct, resStruct)
%    
% -------------------------------------------------------------------------
%  -JJ 1/18/05  change in the format of m3d to use prediction method based on
%  shape space
%   JJ 1/29/08 Following options are Not valid anymore.
%   
%  - writePGA(fileName, header, PGAflag, type, globalStruct, resStruct)
%    => PGAflag = 0 and nFigs > 1   Prediction method based on similarity transform
%  - writePGA(fileName, header, PGAflag, type, globalStruct, resStruct, predStruct)
%    => PGAflag = 0 and nFigs > 1   Prediction method based on shape space
% --------------------------------------------------------------------------------------------
function writePGA(fileName, header, PGAflag, type, globalStruct, resStruct, varargin)

nFigs = sum(header.nFigs);

if (PGAflag == 0)
    numFigStats = 0;
elseif (PGAflag == 1)
    numFigStats = nFigs;
end

isPredOn = 0;
if ( nargin > 6)
    predStruct = varargin{3};
    isPredOn = 1;
end
disp(' ');

nCols = 10;

% ************************************************
% Generate the residue pga statistic file.       *
% ************************************************
fid = fopen(fileName, 'at');

if (fid == -1)
    error(['    Error: cannot open file ' fileName '.']);
end

% ***********************************************
% Generate header                               *
% ***********************************************
fprintf(fid, 'PGAStats {\n');
if (~ismember(type, 0:1))
    type = 1; % mean(r) scaled
end
fprintf(fid, '   scaled = %d;\n', type);

if (isfield(header, 'varScale'))       % CROSS PATIENTS
    fprintf(fid, '   varianceScale = %f;\n', header.varScale);
end

fprintf(fid, '   PGSets {\n');
fprintf(fid, '      numStats = %d;\n', 1+numFigStats);
fprintf(fid, '      Set[0] {\n');
fprintf(fid, '         name = ensemble statistics;\n');
fprintf(fid, '         numFigs = %d;\n', sum(nFigs));
for j = 0:sum(nFigs)-1
    fprintf(fid, '         figIndex[%d] = %d;\n', j, j);
end
fprintf(fid, '      }\n');  % Set[0]

if (PGAflag == 0 || PGAflag == 1)
    for jj = 1:numFigStats
        j = header.dependency.figureOrder(jj);
        fprintf(fid, '      Set[%d] {\n', jj);
        fprintf(fid, ['         name = ' header.figNames{j} ' object residue statistics;\n']);
        fprintf(fid, '         numFigs = %d;\n', 1); % assume single figure objects
        fprintf(fid, '         figIndex[%d] = %d;\n', 0, j-1);
        
        numAugs = 0;
        for k = 1:sum(nFigs)
            if (length(header.dependency.depAtomMatrix{j,k}) > 0)
                numAugs = numAugs + 1;
            end
        end

        fprintf(fid, '         numAugs = %d;\n', numAugs);
        augCnt = 0;
        for k = 1:sum(nFigs)
            numPrims = length(header.dependency.depAtomMatrix{j,k});
            if (numPrims > 0)
                fprintf(fid, '         AugmentedAtoms[%d] {\n', augCnt);
                fprintf(fid, '            figIndex = %d;\n', k-1);
                fprintf(fid, '            numPrims = %d;\n', numPrims);
                for kk = 1:numPrims
                    fprintf(fid, '            primIndex[%d] = %d;\n', kk-1, ...
                        header.dependency.depAtomMatrix{j,k}(kk) - header.figIDs(k));
                end
                fprintf(fid, '         }\n');
                augCnt = augCnt + 1;
            end
        end

        fprintf(fid, '      }\n');  % Set[%d]
    end
end

fprintf(fid, '   }\n'); % PGSets


% ***********************************************
% Ensemble                                      *
% ***********************************************
fprintf(fid, '   Set[0] {\n');

pga = globalStruct.stats.PCs(:, 1:globalStruct.projDims);
PGLength = size(pga,1);

% if (header.adaptive)
%     logMean = atomLog(globalStruct.stats.Mean, 1);
%     logMean = logMean(:);
%     fprintf(fid, '      mean = {%d %d\n', length(logMean), 2);
%     nRows = floor(length(logMean)/nCols);
%     for i = 0:nCols:(nRows-1)*nCols
%         fprintf(fid, '      ');
%         for ii = 1:nCols
%             fprintf(fid, '%12.8f ', logMean(i+ii));
%         end
%         fprintf(fid, '\n');
%     end
%     if (nRows*nCols < length(logMean))
%         fprintf(fid, '      ');
%         for i = nRows*nCols+1:length(logMean)
%             fprintf(fid, '%12.8f ', logMean(i));
%         end
%         fprintf(fid, '\n');
%     end
%     fprintf(fid, '      };\n');
% end

fprintf(fid, '      numPGs = %d;\n', globalStruct.projDims);
fprintf(fid, '      PGLength = %d;\n', PGLength);
for k = 1:globalStruct.projDims
    fprintf(fid, '      PG[%d] = {%d %d \n', k-1, PGLength, 2); % 2 for double
    nRows = floor(PGLength/nCols);
    for i = 0:nCols:(nRows-1)*nCols
        fprintf(fid, '      ');
        for ii = 1:nCols
            fprintf(fid, '%12.8f ', pga(i+ii, k));
        end
        fprintf(fid, '\n');
    end
    if (nRows*nCols < PGLength)
        fprintf(fid, '      ');
        for i = nRows*nCols+1:PGLength
            fprintf(fid, '%12.8f ', pga(i, k));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '      };\n');
end
fprintf(fid, '   }\n');


% ***********************************************
% per Object                                    *
% ***********************************************
 
if (PGAflag == 1)

    for j = 1:numFigStats  

% Move to adapative case
        
%         logMean = mrepLog(resStruct{j}.stats.Mean);
%         [logMeanDataMtx, atomType] = makeDataMtx(logMean);
%         %---------------------
%         % PRINT mean = {}
%         %---------------------
         fprintf(fid, '   Set[%d] {\n', j);
%         fprintf(fid, '      mean = {%d %d\n', length(logMeanDataMtx), 2);
%         nRows = floor(length(logMeanDataMtx)/nCols);
%         for i = 0:nCols:(nRows-1)*nCols
%             fprintf(fid, '      ');
%             for ii = 1:nCols
%                 fprintf(fid, '%12.8f ', logMeanDataMtx(i+ii));
%             end
%             fprintf(fid, '\n');
%         end
%         if (nRows*nCols < length(logMeanDataMtx))
%             fprintf(fid, '      ');
%             for i = nRows*nCols+1:length(logMeanDataMtx)
%                 fprintf(fid, '%12.8f ', logMeanDataMtx(i));
%             end
%             fprintf(fid, '\n');
%         end
%         fprintf(fid, '      };\n');

        % Find out corresponding features of atoms
        pga = makeDataMtx(resStruct(j).stats.PCAtoms);
        pga = pga(1:resStruct(j).projDims, :)';

        PGLength = size(pga, 1);

        %---------------------
        % PRINT PG[%d] = {};
        %---------------------
        fprintf(fid, '      numPGs = %d;\n', resStruct(j).projDims);
        fprintf(fid, '      PGLength = %d;\n', PGLength);
        for k = 1:resStruct(j).projDims
            fprintf(fid, '      PG[%d] = {%d %d\n', k-1, PGLength, 2);
            nRows = floor(PGLength/nCols);
            for i = 0:nCols:(nRows-1)*nCols
                fprintf(fid, '      ');
                for ii = 1:nCols
                    fprintf(fid, '%12.8f ', pga(i+ii, k));
                end
                fprintf(fid, '\n');
            end
            if (nRows*nCols < PGLength)
                fprintf(fid, '      ');
                for i = nRows*nCols+1:PGLength
                    fprintf(fid, '%12.8f ', pga(i, k));
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '      };\n');          % PG[%d] = {};
        end
        %--------------------------------
        % PRINT Prediction{ }
        % JJ 1/30/08: NEED FIX
        %--------------------------------
        if (isPredOn && j < numFigStats)
            fprintf(fid, '     Prediction {\n');
            remAtomIDs = [];
            for jj = j+1:numFigStats
                k = header.dependency.figureOrder(jj);
                remAtomIDs = cat(2, remAtomIDs,  header.figIDs(k)+[0:header.nFigAtoms(k)-1]);
            end
            remLogMean = atomLog(predStruct(j).stats.Mean(:, remAtomIDs), 1);
            remLogMean = remLogMean(:);
            %---------------------
            % PRINT mean = {}
            %---------------------
            fprintf(fid, '          mean = { %d %d\n', length(remLogMean) , 2);
            nRows = floor(length(remLogMean)/nCols);
            for i = 0:nCols:(nRows-1)*nCols
                fprintf(fid, '          ');
                for ii = 1:nCols
                    fprintf(fid, '%12.8f ', remLogMean(i+ii));
                end
                fprintf(fid, '\n');
            end
            if (nRows*nCols < length(remLogMean))
                fprintf(fid, '          ');
                for i = nRows*nCols+1:length(remLogMean)
                    fprintf(fid, '%12.8f ', remLogMean(i));
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '          };\n');  %mean = {

            remPga = predStruct(j).stats.PCs(:,remAtomIDs,1:predStruct(j).projDims);
            remPGLength = size(remPga, 1)*size(remPga, 2);  % dim x nAToms
            remPga = reshape(remPga, remPGLenheader.dependency.figureOrder(jj), predStruct(j).projDims);
            %---------------------
            % PRINT PG[%d] = {};
            %---------------------
            fprintf(fid, '          numPGs = %d;\n', predStruct(j).projDims);
            fprintf(fid, '          PGLength = %d;\n', remPGLength);
            for k = 1:predStruct(j).projDims
                fprintf(fid, '          PG[%d] = {%d %d\n', k-1, remPGLength, 2);
                nRows = floor(remPGLength/nCols);
                for i = 0:nCols:(nRows-1)*nCols
                    fprintf(fid, '          ');
                    for ii = 1:nCols
                        fprintf(fid, '%12.8f ', remPga(i+ii, k));
                    end
                    fprintf(fid, '\n');
                end
                if (nRows*nCols < remPGLength)
                    fprintf(fid, '          ');
                    for i = nRows*nCols+1:remPGLength
                        fprintf(fid, '%12.8f ', remPga(i, k));
                    end
                    fprintf(fid, '\n');
                end
                fprintf(fid, '          };\n');          % PG[%d] = {};
            end

            fprintf(fid, '     }\n');  %Prediction {
        end
        %-----------------------------------
        % END PRINT Prediction{ }
        %-----------------------------------
        fprintf(fid, '   }\n');   %Set[numFigStats]{
    end
end


fprintf(fid, '}\n');    %PGAStats {
fclose(fid);



