function atomList = getAtomListFromSRepData( SRepData, CPNS_Info )     
    
    nRows = CPNS_Info.nRows ;
    nCols = CPNS_Info.nCols ;
    
    nTotalAtoms = nRows * nCols ;
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ;    
    nStdAtoms = nTotalAtoms - nEndAtoms ;    
    
    nTotalPositions = 3 * ( nEndAtoms + nStdAtoms ) ;    
    nTotalRadii = 3 * nEndAtoms + 2 * nStdAtoms ;       
    
    atomList = cell( nRows, nCols ) ;
    
    pVector = SRepData( 1 : nTotalPositions ) ;
    rVector = SRepData( nTotalPositions+1 : nTotalPositions+nTotalRadii ) ;
    uVector = SRepData( nTotalPositions+nTotalRadii+1 : nTotalPositions+nTotalRadii+3*nTotalRadii ) ;   
    
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

            % Get positions

            pIndex = 3 * ( nCols * (i-1) + j-1 ) ;
            
            pos = pVector( pIndex+1 : pIndex+3 ) ;  

            % Get radii

            rIndex = 3 * nEndAtomsBefore + 2 * nStdAtomsBefore ;
            
            if( atomType == 1 )     
                % CASE: END ATOM                 
                r = rVector( rIndex+1 : rIndex+3 ) ;         
            else
                % CASE: STD ATOM
                r = rVector( rIndex+1 : rIndex+2 ) ;                
            end

            % Get spoke directions

            uIndex = 9 * nEndAtomsBefore + 6 * nStdAtomsBefore ;                           

            if( atomType == 1 )                
                % CASE: END ATOM                    
                U = zeros(3, 3) ;
                U( :, 1 ) = uVector( uIndex+1 : uIndex+3 ) ;
                U( :, 2 ) = uVector( uIndex+4 : uIndex+6 ) ;
                U( :, 3 ) = uVector( uIndex+7 : uIndex+9 ) ;                
            else
                % CASE: STD ATOM
                U = zeros(3, 2) ;
                U( :, 1 ) = uVector( uIndex+1 : uIndex+3 ) ;
                U( :, 2 ) = uVector( uIndex+4 : uIndex+6 ) ;
            end                

            % update the no. of atoms

            if( atomType == 1 )
                nEndAtomsBefore = nEndAtomsBefore + 1 ;
            else
                nStdAtomsBefore = nStdAtomsBefore + 1 ;
            end
            
            atomList{i, j} = SrepQuadPrimitive( pos, r, U ) ;
        end
    end    

end