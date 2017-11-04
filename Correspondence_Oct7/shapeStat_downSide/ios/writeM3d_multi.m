
function writeM3d_multi(dataStruct, header, modelType, alignedDataStruct)      

    PROJ.alignedDir='/NIRAL/work/ltu/pablo_matlab/shapeStat/';
            % Number of samples
            nSamps = header.nSamps;
           
            for i = 1:nSamps                
                [pathstr, filename, ext, ver] = fileparts(dataStruct.rawFileNames{i});
                if (isempty(alignedDataStruct.maxDim))
                    writeM3d([PROJ.alignedDir filesep filename '.aligned' ext], ...
                        header, modelType, squeeze(alignedDataStruct.atoms(i, :)));
                else
                    writeM3d([PROJ.alignedDir filesep filename '.aligned' ext], ...
                        header, modelType, squeeze(alignedDataStruct.atoms(i, :)), alignedDataStruct.worldExtents{i});
                end             
               
            end
           
         