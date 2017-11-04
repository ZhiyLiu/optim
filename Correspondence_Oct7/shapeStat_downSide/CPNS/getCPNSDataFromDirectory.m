function [ data info ] = getCPNSDataFromDirectory( dataDir ) 

    %% Read all the models in the data folder %%

    [info m3dData] = readM3d( dataDir, '*.m3d' )  ;
    
    % convert this data to s-reps if they are in old format
    
    if( strcmp( get( m3dData.atoms{1,1}, 'Type' ), 'lie-group' ) )
        disp('CPNS cannot handle lie-group representation of atoms.') ;
        disp('Please convert the models to s-rep representation using Pablo!') ;
        data = [] ;
        return ;
    end    
    
    
    %% Construct the CPNS data matrix by vertical concatenation of p's, r's and U's %%
    
    nSamples = info.nSamps ; 
    
    nRows = info.nRows ;% for me its 3
    nCols = info.nCols ;% its 13
    
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ;    
    nStdAtoms = nRows * nCols - nEndAtoms ;
    
    nTotalPositions = 3 * ( nEndAtoms + nStdAtoms ) ;% each atom has x, y, z.    
    nTotalRadii = 3 * nEndAtoms + 2 * nStdAtoms ;% each spoke has a radii, 106 spokes.    
    nTotalSpokeDirs = 9 * nEndAtoms + 6 * nStdAtoms ;
    
    pMatrix = zeros( nTotalPositions, nSamples ) ;
    rMatrix = zeros( nTotalRadii, nSamples ) ;
    uMatrix = zeros( nTotalSpokeDirs, nSamples ) ;    
    
    for ns = 1 : nSamples,        
        
        nEndAtomsBefore = 0 ;
        nStdAtomsBefore = 0 ;
        
        for i = 1 : nRows,            
            for j = 1 : nCols,
                
                % atomType = 0 ==> std atom
                % atomType = 1 ==> end atom
                
                if( ( i == 1 ) || ( i == nRows ) || ( j == 1 ) || ( j == nCols ) )
                    atomType = 1 ;
                else
                    atomType = 0 ;
                end
                
                atomIndex = nCols * (i-1) + (j-1) + 1 ;
                
                % Get positions
                
                pIndex = 3 * ( nCols * (i-1) + j-1 ) ;
                pMatrix( pIndex+1 : pIndex+3, ns ) = get( m3dData.atoms{ns, atomIndex}, 'pos' ) ;
                
                % Get radii
                
                rIndex = 3 * nEndAtomsBefore + 2 * nStdAtomsBefore ;
                
                rVals = get( m3dData.atoms{ns, atomIndex}, 'r' ) ;
                
                if( atomType == 1 )
                    % CASE: END ATOM                    
                    rMatrix( rIndex+1 : rIndex+3, ns) = rVals' ;
                else
                    % CASE: STD ATOM
                    rMatrix( rIndex+1 : rIndex+2, ns) = rVals(1:2)' ;
                end
                
                % Get spoke directions
                
                uIndex = 9 * nEndAtomsBefore + 6 * nStdAtomsBefore ;                
                
                uVals = get( m3dData.atoms{ns, atomIndex}, 'U' ) ;
                
                if( atomType == 1 )
                    % CASE: END ATOM                    
                    uMatrix( uIndex+1 : uIndex+3, ns) = uVals( :, 1 )' ; 
                    uMatrix( uIndex+4 : uIndex+6, ns) = uVals( :, 2 )' ;
                    uMatrix( uIndex+7 : uIndex+9, ns) = uVals( :, 3 )' ;
                else
                    % CASE: STD ATOM
                    uMatrix( uIndex+1 : uIndex+3, ns) = uVals( :, 1 )' ; 
                    uMatrix( uIndex+4 : uIndex+6, ns) = uVals( :, 2 )' ;
                end                
                
                % update the no. of atoms
                
                if( atomType == 1 )
                    nEndAtomsBefore = nEndAtomsBefore + 1 ;
                else
                    nStdAtomsBefore = nStdAtomsBefore + 1 ;
                end
                
            end
        end
    end
    
    data = [    pMatrix ;
                rMatrix ;
                uMatrix ] ;
end