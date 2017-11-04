% ======================================================================
% Create m3d model file from list of atoms
%   mStruct.atoms is a list of objects.  objID specifies which object
%   to output.  Default objID = 1, and mStruct.atoms has one object.
% ----------------------------------------------------------------------
%modified by xiaoxiao March 2006
%write atomPGA into the m3d file
function writePrimitivePGA(m3dVec,header,atomStruct,fileName, atomNum)

% mStruct.atoms
nFigs = header.nFigs;
nCols = header.nCols;
nRows = header.nRows;
figIDs = header.figIDs;   %where to get the figId????xiaoxiao
nFigAtoms = header.nFigAtoms;


%xiaoxiao :To do: change the format into the same framework of figure PGA.

figTry=1;% atomNum is ordered sequentially according to the figure ID

while( figTry <= size(nFigAtoms,2))
    temp=nFigAtoms(figTry);
    if    ( atomNum <= temp )
        figNo = figTry;
        break;
    else
        atomNum = atomNum-temp;
        figTry = figTry +1;    
    end
end

fid=fopen(fileName,'at');
fprintf(fid,'          deltaAtom[%d] {\n',atomNum-1); 
fprintf(fid,'               atomIndex = %d; \n', atomNum-1);

%mean primitive txt: in order to be used in pablo
% since the residue mean for atom PGA is approximately zero, we don't
% record it for right now
j=figNo;

%gpa statistic

fprintf(fid,'               deltaAtomPGAStats {\n');


%----------------------------------------------------------------------
% Global
%----------------------------------------------------------------------

    %xiaoxiao: pcs is n1Dims*nComp
    pga = atomStruct.stats.PCs(:,1:atomStruct.projDims);
    PGLength = size(pga,1);
    pga = reshape(pga, PGLength, atomStruct.projDims);
    fprintf(fid, '                    numPGs = %d; \n', atomStruct.projDims);
    fprintf(fid, '                    PGLength = %d; \n', PGLength);
    
    nCol=nCols(j);
    
    for k = 1:atomStruct.projDims
        fprintf(fid, '                    PG[%d] = {%d %d \n', k-1, PGLength, 2); % 2 for double
        nRow = floor(PGLength/nCol);
        for i = 0:nCol:(nRow-1)*nCol
            fprintf(fid, '      ');
            for ii = 1:nCol
                fprintf(fid, '%12.8f ', pga(i+ii, k));
            end
          %  fprintf(fid, '\n');
        end
        if (nRow*nCol < PGLength)
            fprintf(fid, '      ');
            for i = nRow*nCol+1:PGLength
                fprintf(fid, '%12.8f ', pga(i, k));
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '                    };\n');
    end

    fprintf(fid, '               }\n');
    fprintf(fid, '          }\n');
    fclose(fid);



