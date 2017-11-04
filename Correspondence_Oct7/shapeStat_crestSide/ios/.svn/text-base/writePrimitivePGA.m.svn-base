% ======================================================================
% Create m3d model file from list of atoms
%   mStruct.atoms is a list of objects.  objID specifies which object
%   to output.  Default objID = 1, and mStruct.atoms has one object.
% ----------------------------------------------------------------------
%modified by xiaoxiao March 2006
%write atomPGA into the m3d file
function writePrimitivePGA(fileName, header,aMatrix,atomNum,atomStruct)

% mStruct.atoms
nFigs = header.nFigs;
nCols = header.nCols;
nRows = header.nRows;
figIDs = header.figIDs;   %where to get the figId????xiaoxiao
nFigAtoms = header.nFigAtoms;
nDims = size(aMatrix, 1);


m3dVec = aMatrix{1};




%xiaoxiao :figId
figTry=1;

while( figTry <= size(nFigAtoms,2))
    temp=nFigAtoms(figTry);
    if    ( atomNum <= temp )
        figNo=figTry;
        figTry=1000;
    else
        atomNum=atomNum-temp;
        figTry = figTry +1;    
    end
end

fid=fopen(fileName,'at');
%fid=fopen([fileName(1:(size(fileName,2)-3)),'txt'],'at');
%mean primitive txt: in order to be used in pablo
fprintf(fid,'          deltaAtom[%d] {\n',atomNum-1); 
fprintf(fid,'               atomIndex = %d; \n', atomNum-1);
fprintf(fid,'               deltaMeanAtom {\n');



j=figNo;
r = floor((atomNum-1) / nCols(j));
c = rem(atomNum-1, nCols(j));
            
            
            if (r==0 | r==nRows(j)-1 | c==0 | c==nCols(j)-1)
                fprintf(fid,'                  elongation = %2.6f;\n',  get(m3dVec,'elongation'));
            end
            


q = get(m3dVec,'q');
fprintf(fid,'                  qw = %2.6f;\n',q(1));
fprintf(fid,'                  qx = %2.6f;\n',q(2));
fprintf(fid,'                  qy = %2.6f;\n',q(3));
fprintf(fid,'                  qz = %2.6f;\n',q(4));


fprintf(fid,'                  r =  %2.6f;\n',get(m3dVec,'r'));
fprintf(fid,'                  theta = %2.6f;\n',get(m3dVec,'theta')*180/pi);

             if (r==0 | r==nRows(j)-1 | c==0 | c==nCols(j)-1)
                fprintf(fid, '                  type = EndPrimitive;\n');
            else
                fprintf(fid, '                  type = StandardPrimitive;\n');
            end
                
p = get(m3dVec,'pos');
fprintf(fid,'                  x = %2.6f;\n',p(1));
fprintf(fid,'                  y = %2.6f;\n',p(2));
fprintf(fid,'                  z = %2.6f;\n',p(3));
fprintf(fid,'                  } \n');

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



