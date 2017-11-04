%===============================================================
% Output eigenvalues at global and residue stage
%==============================================================
function  writeEVs(fileName, header, globalStruct, varargin);

isResidue = 0;
if (nargin == 4)
%    resStruct = varagin;
    resStruct = varargin{1};
%    varargin{1};
    isResidue = 1;
end

fid = fopen(fileName, 'wt');

if (fid == -1)
    disp(['    Error: cannot open file ' fileName '.']);
else
    
    % EVs at Global 
    for i=1:length(globalStruct.stats.EVs)
        fprintf(fid,'%12.8f ', globalStruct.stats.EVs(i));
    end
    fprintf(fid, '\n');
    
    if(isResidue)
        % EVs at residue
        nFigStats = length(header.dependency.figureOrder);
        for i=1:nFigStats
            for j=1:length(resStruct(i).stats.EVs)
                fprintf(fid, '%12.8f ', resStruct(i).stats.EVs(j));
            end
            fprintf(fid, '\n');
        end
    end
end

fclose(fid);