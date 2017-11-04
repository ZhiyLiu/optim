function writeCPNSStats( outputFilename, CPNSShapeStats, CPNSEigenStats, PNSSpoke, scaleSpoke, CPNS_info, MAX_CPNS_MODES ) 

    nRows = CPNS_info.nRows ;
    nCols = CPNS_info.nCols ;
    
    nTotalAtoms = nRows * nCols ;
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ;    
    nStdAtoms = nTotalAtoms - nEndAtoms ;       
    
    nSpokes = 3 * nEndAtoms + 2 * nStdAtoms ;       

    scaleShape  = CPNSShapeStats.scaleShape ;
    meanShape   = CPNSShapeStats.meanShape ;
    PNSShape    = CPNSShapeStats.PNSShape ;   
    
    eigenValues = CPNSEigenStats.eigenValues ;
    eigenVectors = CPNSEigenStats.eigenVectors ;

    fid = fopen( outputFilename, 'a+' ) ;
    
    if( fid < 0 )
        disp('Wrong output filename to writeCPNSStats()') ;
        return ;
    end
    
    % ---------------------------------------------------------------------
    %                  DESCRIPTION OF CPNS STRUCTURE
    
    % PABLO REGISTRY DATA FORMATS
    %
    %     Regular value key: 	keyName = value ;
    %     Array key:            keyName = { ArrayLength ArrayType ArrayValues } ;
    %     Folder key:           keyName { ... }
    %
    %     Array Types = { IntArray (0), FloatArray (1), DoubleArray (2) }
    %
    %     This is how the addition to the CPNS mean file looks like:
    %
    %         CPNS Stats {
    %
    %               nAtomRows = na ;
    %               nAtomCols = nc ;
    % 
    %              scaleShape = gamma ;
    % 
    %              meanShape = { 3,  2, mx, my, mz } ;
    % 
    %              PNSShape {
    %                    nDims  = d ;
    %                     sphereAxis[0] = { d, 2, … } ;
    %                        … 
    %                      sphereAxis[d-1] = vd ;
    %                     sphereDist = { d-1, 2, r1 … rd-1 } ;
    %              } 
    %              nSpokes = ns ;
    % 
    %              scaleSpokes = { ns, 2, rs1 … rsn } ;
    % 
    %              PNSSpoke[0] {
    %              } 
    %              …
    %              PNSSpoke[ns-1] {
    %              }      
    % 
    %              nEigenmodes = ne ;
    % 
    %              eigenValues = { ne, 2, ... }  ;
    % 
    %              eigenVector[0] = { lv, 2, … } ;
    %                    …
    %              eigenVector[ne] = { lv, 2, … } ;
    % 
    %         }  # end of CPNSStats
    
    
    % ---------------------------------------------------------------------
    
    fprintf( fid, 'CPNSStats { \n' ) ;
    
        fprintf( fid, 'nAtomRows = % d ; \n', nRows ) ;
        
        fprintf( fid, 'nAtomCols = % d ; \n', nCols ) ;
    
        fprintf( fid, 'scaleShape = % 12.8f ; \n', scaleShape ) ;

        fprintf( fid,   'meanShape = { %d %d % 12.8f % 12.8f % 12.8f } ; \n', ... 
                        3, 2, meanShape(1), meanShape(2), meanShape(3) ) ;

        % ----------- PNS Shape -------------------------------------------

        fprintf( fid, 'PNSShape { \n' ) ;

            nDims = length( PNSShape.radii ) ;

            fprintf( fid, 'nDims = %d ; \n', nDims ) ;

            for d = 1 : nDims,
                
                if( d < nDims )
                    fprintf( fid, 'sphereAxis[%d] = { %d %d', d-1, nDims-(d-1)+1, 2 ) ;
                    for i = 1 : length( PNSShape.orthaxis{d} ),
                        fprintf( fid, ' %12.8f ', PNSShape.orthaxis{d}(i) ) ;
                    end
                    fprintf( fid, '} ; \n' ) ;  % end of sphereAxis
                else
                    fprintf( fid, 'sphereAxis[%d] = { %d %d', d-1, 1, 2 ) ;                        
                    fprintf( fid, ' %12.8f ', PNSShape.orthaxis{d}(1) ) ;                        
                    fprintf( fid, '} ; \n' ) ;  % end of sphereAxis                                            
                end
            end

            fprintf( fid, 'sphereDist = { %d %d ', nDims-1, 2 ) ;

            for d = 1 : nDims-1,
                fprintf( fid, ' %12.8f ', PNSShape.dist(d) ) ;
            end

            fprintf( fid, '} ; \n' ) ;  % end of sphereDist

        fprintf( fid, '} \n' ) ;  % end of PNSShape
        
        % -------------------------------------------------------
        
        fprintf( fid, 'nSpokes = %d ; \n', nSpokes ) ;
        
        fprintf( fid, 'scaleSpoke = { %d %d ', nSpokes, 2 ) ;

        for n = 1 : nSpokes,
            fprintf( fid, ' %12.8f ', scaleSpoke(n) ) ;
        end        
        
        fprintf( fid, '} ; \n' ) ;  % end of scaleSpoke
        
        % ----------- PNS Spokes ------------------------------------------
        
        for n = 1 : nSpokes,

            fprintf( fid, 'PNSSpoke[%d] { \n', n-1 ) ;

                nDims = length( PNSSpoke{n}.radii ) ;

                fprintf( fid, 'nDims = %d ; \n', nDims ) ;

                for d = 1 : nDims,                    
                    if( d < nDims )                        
                        fprintf( fid, 'sphereAxis[%d] = { %d %d', d-1, nDims-(d-1)+1, 2 ) ;
                        for i = 1 : length( PNSSpoke{n}.orthaxis{d} ),
                            fprintf( fid, ' %12.8f ', PNSSpoke{n}.orthaxis{d}(i) ) ;
                        end
                        fprintf( fid, '} ; \n' ) ;  % end of sphereAxis
                    else
                        fprintf( fid, 'sphereAxis[%d] = { %d %d', d-1, 1, 2 ) ;                        
                        fprintf( fid, ' %12.8f ', PNSSpoke{n}.orthaxis{d}(1) ) ;                        
                        fprintf( fid, '} ; \n' ) ;  % end of sphereAxis                        
                    end
                end   

                fprintf( fid, 'sphereDist = { %d %d ', nDims-1, 2 ) ;

                for d = 1 : nDims-1,
                    fprintf( fid, ' %12.8f ', PNSSpoke{n}.dist(d) ) ;
                end

                fprintf( fid, '} ; \n' ) ;  % end of sphereDist

            fprintf( fid, '} \n' ) ;  % end of PNSSpoke{n}
            
        end        
        
        % ----------------- EigenVectors and EigenValues ------------------
        
        fprintf( fid, 'nEigenmodes = %d ; \n', MAX_CPNS_MODES ) ;
        
        fprintf( fid, 'eigenValue = { %d %d ', MAX_CPNS_MODES, 2 ) ;

        for n = 1 : MAX_CPNS_MODES,
            fprintf( fid, ' %12.8f ', eigenValues(n) ) ;
        end

        fprintf( fid, '} ; \n' ) ;  % end of eigenValues
        
        l = size(eigenVectors, 1) ; % length of each eigenVector
        
        fprintf( fid, 'eigenVectorLength = %d ; \n', l ) ;        
        
        for n = 1 : MAX_CPNS_MODES,
            fprintf( fid, 'eigenVector[%d] = { %d %d ', n-1, l, 2 ) ;

            for i = 1 : l,
                fprintf( fid, ' %12.8f ', eigenVectors( i, n) ) ;
            end
            fprintf( fid, '} ; \n' ) ;  % end of eigenVector
        end           
    
        % -----------------------------------------------------------------
        
    fprintf( fid, '} \n' ) ;    % end of CPNSStats
    
    fclose( fid ) ;    

end