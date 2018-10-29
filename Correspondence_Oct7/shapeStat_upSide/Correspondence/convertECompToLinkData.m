function linkData = convertECompToLinkData(CPNS_info, CPNSScores, meanSizeOfPDMs, meanOfCombinedPDM, meanRs, PNSShape, PNSSpoke )
    nLinks = CPNS_info.nLinks;
    nSamples = CPNS_info.nSamples;
    %% Construct the S-rep data (X, radii, spokeDirs) from Zcomp (CPNSScores)

    % 1. Construct X (hub positions) 
    %%%%%%% check for dimension > sample size case (3/25/2013 by S Jung)
%     cpnsdim = size(CPNSScores,1);
%     if cpnsdim < nTotalPositions + 3 * nTotalRadii; 
%        dimZShape = cpnsdim - 3 * nTotalRadii;
%     else
%        dimZShape = nTotalPositions;
%     end
    %%%%%%%
    
    dimZShape = nLinks;
    
    zShape = CPNSScores( 1 : dimZShape-1 ) / meanSizeOfPDMs;     % = zShape part of z vector, undo pre-PCA scaling

    zSizePDM = CPNSScores( dimZShape );   % = z_gamma

    sizePDMOverall = meanSizeOfPDMs * exp( zSizePDM / meanSizeOfPDMs );     % = gamma
    
    XStar = PNSe2s( zShape, PNSShape ) ;  % = X*
        
    X = sizePDMOverall * XStar + repmat( meanOfCombinedPDM', [ nTotalAtoms, 1] ) ;  % = gamma x X* + X_bar
   
end