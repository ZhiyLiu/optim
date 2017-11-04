function calculateCPNS( dataDir, outputFilename, optionPNS )
% calculateCPNS:    Function to calculate the CPNS mean and modes for a
%                   collection of s-rep models in a folder. This program
%                   writes out the mean model and the eigenmode information
%                   in a "*.m3d" file which can be read by Binary Pablo.
%                   The program also outputs a "*.jpeg" image containing 
%                   the plot of the eigenvalues (first 20)
%
% INPUT PARAMETERS: calculateCPNS( dataDir, outputFilename, optionPNS )
%
% dataDir => The folder where the s-reps are stored. This is also the
%               folder where the output mean file and the eigenvalues plot
%               is written
%
% outputFilename => The name of the output mean file. If this includes a
%                   folder, then the output mean file is written to that 
%                   folder. Otherwise it is written to the dataDir.  
% 
% optionPNS => specify how CPNS should operate
%
% optionPNS = 0 => use SMALL circles or BIG circles depending upon analysis
%           = 1 => always use SMALL circles
%           = 2 => always use BIG circles
%
% CONTROL FLAGS (to be set inside the program code)
%
%   MAX_CPNS_MODES          => maximum no. of modes to be written (~10-15)
%   subtractNoisyVariance   => whether to subtract noisy variance {0,1}
%   verbosity               => whether comments are displayed or not {0,1}
%   loadCPNSModelData       => whether to load previously saved CPNS model data {0,1}
%   saveCPNSModelData       => whether to save model data after loading {0,1}
%   loadStoredPNSShape      => whether to load previously saved PNSShape {0,1}
%   loadStoredPNSSpokes     => whether to load previously saved PNSSpokes {0,1}
%   saveComputedPNS         => whetehr to save computed PNSShape and PNSSpokes {0,1}


if( nargin < 3 )
    optionPNS = 0 ;
end

%%                          CONTROL FLAGS                               %%

% The maximum no. of CPNS modes to be written to the CPNS mean file
MAX_CPNS_MODES = 15 ;   

% subtractNoisyVariance: This flag, when ON, implies that we consider the 
%      variance after the MAX_CPNS_MODES to be noise (because of 
%      irregularities in observations and error in pre-alignment of 
%      samples). So we subtract the average variance (average of 
%      eigenValues) of the remaining modes from the first MAX_CPNS_MODES 
%      no. of eigenModes
subtractNoisyVariance = 0 ;

% whether to display status messages or not
verbosity = 1 ;         

% loadCPNSModelData     =>  specify whether CPNS model data should be
%                           loaded from a previously saved .mat file
% loadCPNSModelData     = 1 =>  load CPNS model data from .mat file
%                       = 0 =>  CPNS data is produced by reading .m3d 
%                               models from the data folder
loadCPNSModelData = 1 ;

% saveCPNSModelData     =>  specify whether CPNS model data (from .m3d files) 
%                           should saved in .mat file to be used later
% saveCPNSModelData     = 1 => save data into the .mat file
saveCPNSModelData = 1 ;

% loadStoredPNS*    =>  specify whether CPNS data should be loaded from a
%                       previously saved .mat file
% loadStoredPNS*    = 1 => load data from .mat file
%                   = 0 => PNS would be computed every time
loadStoredPNSShape = 1 ;
loadStoredPNSSpokes = 1 ;

% saveComputedPNS   =>  specify whether computed CPNS data should saved in 
%                       .mat file to be used later
% saveComputedPNS   = 1 => save data into the .mat file
saveComputedPNS = 1 ;


%%                    Parse input arguments                             %%

    % if no directory is specified in the outputFilename, the output mean 
    % and eigenValues plot are saved in the data directory itself
    
    [pathstr, name, ext] = fileparts(outputFilename) ;
    
    if( isempty( pathstr) )
        outputFilename = fullfile( dataDir, [name '.m3d'] ) ;
    else
        outputFilename = fullfile( pathstr, [name '.m3d'] ) ;
    end
    
%%  Get data from the data folder containing fitted s-rep models %%    

    readSRepModels = 1 ;

    if( loadCPNSModelData )
        CPNSModelDataFile = fullfile( dataDir, 'CPNSModelData.mat' ) ;
        if( exist( CPNSModelDataFile, 'file' ) )
            load( CPNSModelDataFile ) ;
            readSRepModels = 0 ;
            disp(['CPNS s-rep models data loaded from file: ' CPNSModelDataFile]) ;
        end
    end
    
    if( readSRepModels )            
        [ CPNS_data CPNS_info ] = getCPNSDataFromDirectory( dataDir ) ;        
        
        if( isempty(CPNS_data) )
            disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
            return ;
        else        
            disp(['CPNS Data has been read from the folder: ' dataDir]) ;    
        end        
    end
    

    
    if( saveCPNSModelData )
        CPNSModelDataFile = fullfile( dataDir, 'CPNSModelData.mat' ) ;
        save( CPNSModelDataFile, 'CPNS_data', 'CPNS_info' ) ;
        disp(['CPNS s-rep models data saved into file: ' CPNSModelDataFile]) ;
    end
    
    nSamples = CPNS_info.nSamps ; 
    
    nRows = CPNS_info.nRows ;
    nCols = CPNS_info.nCols ;
    
    nTotalAtoms = nRows * nCols ;
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ;    
    nStdAtoms = nTotalAtoms - nEndAtoms ;    
    
    nTotalPositions = 3 * ( nEndAtoms + nStdAtoms ) ;    
    nTotalRadii = 3 * nEndAtoms + 2 * nStdAtoms ;       
       
       
%% CPNS: Step 1 : Deal with hub Positions (PDM) %%    

    % position(i, j, k) =  j-th co-ordinate of i-th atoms position of the kth sample
    position = zeros( nTotalAtoms, 3, nSamples ) ;   

    for i = 1 : nTotalAtoms,
        for j = 1 : 3,                
            position( i, j, : ) = CPNS_data( 3*(i-1)+j, : ) ;
        end
    end     
    
    meanOfEachPDM = mean( position, 1 ) ;                  % 1X3XN matrix storing the mean of PDM of each sample
    
    meanOfCombinedPDM = mean( mean( position, 1 ), 3 ) ;    % 1X3 matrix storing the mean of entire PDM (all samples)
    
    % translated positions centered at the mean of each sample
    
    cposition = position - repmat( meanOfEachPDM, [nTotalAtoms, 1, 1] ) ;
    
    % cposition = position - repmat( meanOfCombinedPDM, [nTotalAtoms, 1, nSamples] ) ; 

    % dibyendu
    % sscarryPDM = ( 1 X 1 X nSamples ) array for holding size of each PDM
    % sscarryPDM(k) = size of the PDM for k-th sample 
    
    sscarryPDM = sqrt(sum(sum(cposition.^2))); 
    
    % sphmatPDM = position part of the CPNS data (translated + scaled)
    % this is the input to the PNS program for the SHAPE part
    
    sphmatPDM = zeros( nTotalPositions, nSamples );
    spharrayPDM = zeros( nTotalAtoms, 3, nSamples );

    for i = 1 : nSamples,
        spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
        sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nTotalAtoms * 3, 1 );
    end

    % ------------------ Shape analysis via PNS ------------------     
    
    % load the stored PNS data from a .mat file if it has already been
    % computed
    
    computePNSShape = 1 ;
    
    % load stored data if it exists
    if( loadStoredPNSShape )        
        PNSShapeDataFile = fullfile( dataDir, 'PNSShape.mat');        
        if ( exist( PNSShapeDataFile, 'file' ) )            
            load( PNSShapeDataFile ) ;
            computePNSShape = 0 ;
            disp(['PNS for shape (PDM) loaded from file: ' PNSShapeDataFile]) ;
        end        
    end
    
    % compute PNS
    
    if( computePNSShape )
        [ZShape PNSShape] = PNSmain( sphmatPDM );  
        %[ZShape PNSShape] = PNSmain( sphmatPDM, optionPNS, 0.05, 100, true);
       
    end
    
    debugMode = 0 ;
    
    if( debugMode )        
        E1 = ZShape(:,1) ;
        S1 = sphmatPDM(:,1) ;        
    end
    
    
    % save PNS data if required    
    if( saveComputedPNS )
        PNSShapeDataFile = fullfile( dataDir, 'PNSShape.mat');
        save( PNSShapeDataFile, 'ZShape', 'PNSShape' ) ;
        disp(['PNS for shape (PDM) saved into file: ' PNSShapeDataFile]) ;
    end
    
    if( verbosity )
        disp('PNS computed for Shape PDMs') ;
    end    
    
    % Size of PDM's 
    
    sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
    sizePDM(:) = sscarryPDM(1, 1, :) ;    

    meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

    normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables

%% CPNS: Step 2 : Deal with atom radii (log radii are Euclidean variables) %%    

    logR = zeros( nTotalRadii, nSamples ) ;     % logR(i, j): i-th radius of j-th sample
    
    for i = 1 : nTotalRadii,        
        logR(i,:) = log(CPNS_data(i+nTotalPositions,:));
    end

    meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
    meanRs = exp( meanLogR ) ;                  % r_bar(i): mean of i-th radius along all samples - stored as scaleSpoke
    
    % r_bar(i, :): mean of i-th radius for all samples
    % used as a scaling factor for radii and spoke directions
    rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
    
    uScaleFactors =  zeros( 2 * nTotalRadii, nSamples ) ;
    for i = 1 : nTotalRadii,
        uScaleFactors( 2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
        uScaleFactors( 2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
    end
   

    
    % shifted log r-i's
    RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)    
    

%% CPNS: Step 3 : Deal with spoke directions (direction analysis) %% %% 

    computePNSSpokes = 1 ;
    
    % load stored data if it exists
    if( loadStoredPNSSpokes )        
        PNSSpokesDataFile = fullfile( dataDir, 'PNSSpokes.mat');        
        if ( exist( PNSSpokesDataFile, 'file' ) )            
            load( PNSSpokesDataFile ) ;
            computePNSSpokes = 0 ;
            disp(['PNS for spokes (direction analysis) loaded from file: ' PNSSpokesDataFile]) ;
        end        
    end

    if( computePNSSpokes )
        
        PNSSpoke = cell( nTotalRadii, 1 ) ;

        ZSpoke = zeros( 2 * nTotalRadii,  nSamples ) ;   

        tOpStart = tic ;

        for r = 1 : nTotalRadii,

            a = CPNS_data( nTotalPositions+nTotalRadii+(r-1)*3+1 : nTotalPositions+nTotalRadii+(r-1)*3+3,:);

            [Zr PNSr] = PNSmain(a,optionPNS);
            %[Zr PNSr] = PNSmain(a,optionPNS, 0.05, 100, true);

%            tElapsed = toc( tOpStart ) ;
%            tRemaining = ( tElapsed / r ) * ( nTotalRadii - r ) ;

%             if( verbosity )
%                 disp( ['Computing PNS for spoke ' num2str(r) ' of ' num2str(nTotalRadii)] ) ;
%                 if( tRemaining > 60 )
%                     disp( ['Estimated time remaining is ' num2str(tRemaining/60) ' mins'] ) ;
%                 else
%                     disp( ['Estimated time remaining is ' num2str(tRemaining) ' secs'] ) ;
%                 end
%                 if( r == nTotalRadii )
%                     disp(['Total time elasped for PNS calculation of ' num2str(nTotalRadii) ' spokes = ' num2str(tElapsed/60) ' mins']) ;
%                 end
%             end                    

            PNSSpoke{r} = PNSr;

            ZSpoke( 2*(r-1) + 1 : 2*(r-1) + 2, : ) = Zr ;
        end        
    end
    
    % save PNS data if required    
    if( saveComputedPNS )
        PNSSpokesDataFile = fullfile( dataDir, 'PNSSpokes.mat');
        save( PNSSpokesDataFile, 'ZSpoke', 'PNSSpoke' ) ;
        disp(['PNS for spokes (direction analysis) saved into file: ' PNSSpokesDataFile]) ;
    end    

%% CPNS: Step 4 : Construct composite Z matrix and perform PCA %% %%

    ZComp = [   meanSizePDM * ZShape ; ... 
                meanSizePDM * normalizedSizePDM ; ... 
                rScaleFactors .* RStar ; ... 
                uScaleFactors .* ZSpoke  ];
            
    [ eigenVectors, pScores, eigenValues ]= princomp( ZComp', 'econ' ) ;
    
    nEigenModesToWrite = min( MAX_CPNS_MODES, length(eigenValues) ) ;
    
    if( subtractNoisyVariance )        
        if( nEigenModesToWrite < length( eigenValues ) )
            eigenValToSubtract = mean( eigenValues(nEigenModesToWrite+1:end) ) ;
            eigenValues = eigenValues - eigenValToSubtract ;
            eigenValues( eigenValues < 0 ) = 0.0 ;
        end        
    end
    

%% Plot the graph for PCA proportions and save as an image %%

    nPlots = min( length( eigenValues ), 2 * MAX_CPNS_MODES ) ;

    hFig = figure(1); clf;
    eigen_prop = eigenValues/sum(eigenValues)*100 ;
    bar(eigen_prop(1:nPlots)) ;
    hold on;
    cum_prop = cumsum(eigenValues)/sum(eigenValues)*100; 
    
    for i=1:length(cum_prop)
        fprintf('%d: %f\n', i, cum_prop(i));
    end
    
    plot( cum_prop(1:nPlots) ) ;
    xlim( [1 nPlots] ) ;
    ylim( [0 100] ) ;
    title('Plot of Eigenmode contributions for CPNS') ;
    xlabel( 'No. of Eigenmodes' ) ;
    ylabel('Percentage contribution') ;   
    
    [pathstr, name, ext] = fileparts(outputFilename) ;
    
    imgEigenValsFilename = fullfile(pathstr, [name '_CPNS_eVals.jpg']) ;    
    saveas( hFig, imgEigenValsFilename ) ;
    disp(['Eigenvalues image written to: ' imgEigenValsFilename ]) ;
    
    
%% CPNS: Step 5: Transform mean back from E-comp (eucliden) space to S (SRep) space %% 
    
    % 1. Convert from EComp to S-space
    
        CPNSMeanScores = zeros( size(eigenVectors, 1), 1 ) ;     % = z vector    
        
        SRepData = convertECompToSRepData( CPNS_info, CPNSMeanScores, meanSizePDM, meanOfCombinedPDM, meanRs, PNSShape, PNSSpoke ) ;     
        
    % 2. Convert SRepData to an s-rep (Quad)
    
        atomList = getAtomListFromSRepData( SRepData, CPNS_info ) ;
    
    % 3. Write out the mean and the CPNS statistics
    
        % writing the mean model
        
        modelType = 1 ;     % mean
    
        writeM3d( outputFilename, CPNS_info, modelType, reshape(atomList', 1, nRows*nCols) ) ;    

        % writing the CPNS Statistics
        
        % packaged into separate struct's because Matlab's rule that all
        % cell arrays in the same struct should be of the same size
        
        CPNSShapeStats = struct(    'scaleShape',   meanSizePDM , ...
                                    'meanShape',    meanOfCombinedPDM, ...
                                    'PNSShape',     PNSShape    ) ;
                                
        CPNSEigenStats = struct(    'eigenValues',  eigenValues, ...
                                    'eigenVectors', eigenVectors        ) ;                       
                        
        scaleSpoke = meanRs ;
                                
        % write CPNS Stats
        
        writeCPNSStats( outputFilename, CPNSShapeStats, CPNSEigenStats, PNSSpoke, scaleSpoke, CPNS_info, nEigenModesToWrite ) ;
        
        close 'all' ;
        
end
