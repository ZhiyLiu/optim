
function SRepData = convertECompToSRepData( CPNS_info, CPNSScores, meanSizeOfPDMs, meanOfCombinedPDM, meanRs, PNSShape, PNSSpoke )

    %% Get info about the S-rep from CPNSInfo %%
    
    nRows = CPNS_info.nRows ;
    nCols = CPNS_info.nCols ;
    
    nTotalAtoms = nRows * nCols ;
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ;    
    nStdAtoms = nTotalAtoms - nEndAtoms ;    
    
    nTotalPositions = 3 * ( nEndAtoms + nStdAtoms ) ;    
    nTotalRadii = 3 * nEndAtoms + 2 * nStdAtoms ;       
    
    %% Construct the S-rep data (X, radii, spokeDirs) from Zcomp (CPNSScores)

    % 1. Construct X (hub positions) 
    
        zShape = CPNSScores( 1 : nTotalPositions-1 );     % = zShape part of z vector

        zSizePDM = CPNSScores( nTotalPositions );   % = z_gamma

        sizePDMOverall = meanSizeOfPDMs * exp( zSizePDM / meanSizeOfPDMs );     % = gamma

        XStar = PNSe2s( zShape, PNSShape ) ;                                    % = X*               

        X = sizePDMOverall * XStar + repmat( meanOfCombinedPDM', [ nTotalAtoms, 1] ) ;  % = gamma x X* + X_bar
    
    % 2. Construct radii (r's)
    
        zRStar = CPNSScores( nTotalPositions + 1 : nTotalPositions + nTotalRadii );     % = z_R* [nSpokes]
        
        radii = meanRs .* exp( zRStar ./ meanRs );
    
    % 3. Construct spokeDirs (U's)
    
        zSpokes = CPNSScores( nTotalPositions + nTotalRadii + 1 : nTotalPositions + nTotalRadii + 2*nTotalRadii );    % = z_Spokes [2 X nSpokes]
    
        spokeDirs = zeros( 3 * nTotalRadii, 1 ) ;

        for ns = 1 : nTotalRadii,
            spokeDirs( 3*(ns-1) + 1 : 3*(ns-1) + 3 ) = PNSe2s( zSpokes(2*(ns-1)+(1:2),:), PNSSpoke{ns} );
        end

        SRepData = [X ; radii ; spokeDirs];

end