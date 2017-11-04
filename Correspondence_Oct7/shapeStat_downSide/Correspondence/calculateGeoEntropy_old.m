% side: 0, up spokes; 1, down spokes.
%function [ geoEntropy] = calculateGeoEntropy( dataDir, outputFilename, optionPNS, side )
function [ geoEntropy] = calculateGeoEntropy( dataDir, optionPNS, side )
% dataDir = '/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/withdeltauv';
% outputFilename='output-plot_LV';
% optionPNS = 1;
% side =0;

if( nargin < 3 )
    optionPNS = 0 ;
end

%%                          CONTROL FLAGS                               %%

% subtractNoisyVariance: This flag, when ON, implies that we consider the 
%      variance after the MAX_CPNS_MODES to be noise (because of 
%      irregularities in observations and error in pre-alignment of 
%      samples). So we subtract the average variance (average of 
%      eigenValues) of the remaining modes from the first MAX_CPNS_MODES 
%      no. of eigenModes

% whether to display status messages or not
verbosity = 1 ;         

% loadCPNSModelData     =>  specify whether CPNS model data should be
%                           loaded from a previously saved .mat file

% saveCPNSModelData     =>  specify whether CPNS model data (from .m3d files) 
%                           should saved in .mat file to be used later
% saveCPNSModelData     = 1 => save data into the .mat file
saveCPNSModelData = 1 ;

% loadStoredPNS*    =>  specify whether CPNS data should be loaded from a
%                       previously saved .mat file
% loadStoredPNS*    = 1 => load data from .mat file
%                   = 0 => PNS would be computed every time
loadStoredPNSShape = 1 ;

% loadCPNSModelData     = 1 =>  load CPNS model data from .mat file
%                       = 0 =>  CPNS data is produced by reading .m3d 
%                               models from the data folder
loadCPNSModelData = 1 ;

% saveComputedPNS   =>  specify whether computed CPNS data should saved in 
%                       .mat file to be used later
% saveComputedPNS   = 1 => save data into the .mat file
saveComputedPNS = 1 ;


%%                    Parse input arguments                             %%

    % if no directory is specified in the outputFilename, the output mean 
    % and eigenValues plot are saved in the data directory itself
    
   % [pathstr, name, ext] = fileparts(outputFilename) ;
    
   % if( isempty( pathstr) )
   %     outputFilename = fullfile( dataDir, [name '.m3d'] ) ;
   % else
   %     outputFilename = fullfile( pathstr, [name '.m3d'] ) ;
   % end
    
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
        [ CPNS_data CPNS_info ] = getCPNSDataMatrix( dataDir, side ) ;        
        
        if( isempty(CPNS_data) )
            disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
            return ;
        else        
            disp(['CPNS Data has been read from the folder: ' dataDir]) ;    
        end        
    end    

    % save CPNS_data and CPNS_info into .mat file
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
    %nTotalRadii = 3 * nEndAtoms + 2 * nStdAtoms ;
    nTotalUpRadii =  nEndAtoms + nStdAtoms;
 
%% CPNS: Step 1 : Deal with hub Positions (PDM) %%    

    % position(i, j, k) =  j-th co-ordinate of i-th atoms position of the kth sample
    position = zeros( nTotalAtoms, 3, nSamples ) ;   

    for i = 1 : nTotalAtoms,
        for j = 1 : 3,                
            position( i, j, : ) = CPNS_data( 3*(i-1)+j, : ) ;
        end
    end     
    
    meanOfEachPDM = mean( position, 1 ) ;% 1X3XN matrix storing the mean of PDM of each sample    
   
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
        [ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100 );  %add default parameters values.        
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

    logR = zeros( nTotalUpRadii, nSamples ) ;     % logR(i, j): i-th radius of j-th sample
    
    for i = 1 : nTotalUpRadii,        
        logR(i,:) = log(CPNS_data(i+nTotalPositions,:));
    end

    meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
    
    % r_bar(i, :): mean of i-th radius for all samples
    % used as a scaling factor for radii and spoke directions
    rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
    
    uScaleFactors =  zeros( 2 * nTotalUpRadii, nSamples ) ;
    for i = 1 : nTotalUpRadii,
        uScaleFactors( 2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
        uScaleFactors( 2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
    end
    
    % shifted log r-i's
    RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  

%% CPNS: Step 3 : Deal with spoke directions (direction analysis) %% %%
       
        PNSSpoke = cell( nTotalUpRadii, 1 ) ;

        ZSpoke = zeros( 2 * nTotalUpRadii,  nSamples ) ; 

        tOpStart = tic ;

        for r = 1 : nTotalUpRadii,

            a = CPNS_data( nTotalPositions+nTotalUpRadii+(r-1)*3+1 : nTotalPositions+nTotalUpRadii+(r-1)*3+3,:);

            [Zr PNSr] = PNSmain(a,optionPNS);

            tElapsed = toc( tOpStart ) ;
            tRemaining = ( tElapsed / r ) * ( nTotalUpRadii - r ) ;

            if( verbosity )
                disp( ['Computing PNS for spoke ' num2str(r) ' of ' num2str(nTotalUpRadii)] ) ;
                if( tRemaining > 60 )
                    disp( ['Estimated time remaining is ' num2str(tRemaining/60) ' mins'] ) ;
                else
                    disp( ['Estimated time remaining is ' num2str(tRemaining) ' secs'] ) ;
                end
                if( r == nTotalUpRadii )
                    disp(['Total time elasped for PNS calculation of ' num2str(nTotalUpRadii) ' spokes = ' num2str(tElapsed/60) ' mins']) ;
                end
            end                    

            PNSSpoke{r} = PNSr;
           
            ZSpoke( 2*(r-1) + 1 : 2*(r-1) + 2, : ) = Zr ;
             %disp(ZSpoke);
        end
        %disp(ZSpoke);
    
    % save PNS data if required    
    if( saveComputedPNS )
        PNSSpokesDataFile = fullfile( dataDir, 'PNSSpokes.mat');
        save( PNSSpokesDataFile, 'ZSpoke', 'PNSSpoke') ;
        disp(['PNS for spokes (direction analysis) saved into file: ' PNSSpokesDataFile]) ;
    end    

%% CPNS: Step 4 : Construct composite Z matrix and perform PCA %% %%

    ZComp = [   meanSizePDM * ZShape ; ... 
                meanSizePDM * normalizedSizePDM ; ... 
                rScaleFactors .* RStar ; ... 
                uScaleFactors .* ZSpoke  ];

%% princomp method use corvariance matrix to calculate the eigenvalue. This
% output the same result as the following standard PCA method.
%    [ eigenVectors, pScores, eigenValues ]= princomp( ZComp', 'econ' ) ;
%    disp('The eigenVectors is:');
%    disp(eigenVectors);
%    
%    disp('The pScores is:');
%     disp(pScores);
%     
%     disp('The eigenvalues is:');
%     disp(eigenValues);
    
%% Standard PCA
% % ZComp is a matrix with row as feature(dimesion), column as a case(sample).
% [m,n] = size(ZComp);
% 
% % mean for each dimension(one row is one dimension), its a m by 1 matrix.
% C = mean(ZComp,2);
% 
% % subtract off the row mean from A
% M = repmat(C,1,n); %repmat(C,1,n) means tile 1*n patch of C.
% centereddata = ZComp - M;
% 
% % calculate the covariance matrix
% trancentereddata = transpose(centereddata);
% covariance_1 = 1/(n-1) * centereddata * trancentereddata;
% % make the covariance matrix symmetric.
% covariance = (covariance_1 + transpose(covariance_1))/2;
% [PC, V] = eig(covariance);
% 
% % extract diagonal of matrix (eigenvalues) as vector
% V = diag(V);
% 
% % sort the variances in decreasing order
% % rindices is the original index of element in V, use this we can get the
% % element from V in decreasing order.
% [junk, rindices] = sort(-1*V);
% eigenvalue = V(rindices);
% disp(eigenvalue);

%% PCA using correlation matrix instead of corvariance matrix
% ZComp is a matrix with row as feature(dimesion), column as a case(sample).
[m,n] = size(ZComp);

% mean for each dimension(one row is one dimension), its a m by 1 matrix.
C = mean(ZComp,2);

% subtract off the row mean from A
M = repmat(C,1,n); %repmat(C,1,n) means tile 1*n patch of C.
centereddata = ZComp - M;

% delete the rows if its elements all zero.
zindex = all(centereddata==0,2);
%disp('The row index is: ');
%disp(zindex);
centereddata(zindex,:) = [];

% calculate the covariance matrix
trancentereddata = transpose(centereddata);

covariance_1 = 1/(n-1) * centereddata * trancentereddata;
% make the covariance matrix symmetric.
covariance = (covariance_1 + transpose(covariance_1))/2;
disp('The covariance matrix for entropy is: ');
disp(covariance);

% Convert covariance matrix to correlation matrix.
correlation_1 = corrcov(covariance);
correlation = (correlation_1 + transpose(correlation_1))/2;

% find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V)
[PC, V] = eig(correlation);

% extract diagonal of matrix (eigenvalues) as vector
V = diag(V);
    
% descent order
[junk, rindices] = sort(-1*V);
eigenvalue = V(rindices);
[eRows, eCols] = size(eigenvalue);
indexi = 1;
for ind = 1:eRows
    if(eigenvalue(ind,1)>1e-6)  %only keep the eigenvalue which is big than 1e-6.
        EigenvalueCut(indexi,1) = eigenvalue(ind,1);
        indexi = indexi +1;
    end
end
%disp('The eigen value cuted is: ');
%disp(EigenvalueCut);

% Calculate the entropy: 
geoEntropy = 0;
for i = 1:(indexi-1)    %only use the eigenvalues big than 1e-6.
    geoEntropy = geoEntropy + log(eigenvalue(i,1));
end  
%  disp('The geometry entropy is: ');
%  disp(geoEntropy);
    
