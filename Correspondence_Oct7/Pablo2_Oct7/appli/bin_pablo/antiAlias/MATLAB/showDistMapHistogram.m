function showDistMapHistogram( raw3File, dataType , imgSize)

    if( nargin < 3 )
        imgSize = 128 * 128 * 128 ;
    end
    
    if( nargin < 2 )
        dataType = 'ushort' ;
    end

    [header data] = readRaw3()    ;
    
    hist( data, 1000 ) ;
    
    display('Perform operations on the data part now') ;
    
    close all ;
    
    function [header1, data1] = readRaw3()

        fid1 = fopen(raw3File, 'r');

        switch( dataType )

            case 'ushort'    
                fileData = fread(fid1, 'uint16');    

            case 'short'
                fileData = fread(fid1, 'int16');    

        end   

        nHeader = size(fileData,1) - imgSize ;

        header1 = fileData(1:nHeader) ;    

        data1 = fileData(nHeader+1:size(fileData,1)) ;    

        fclose(fid1);
    end
    
    
end