function [ ] = geoEntropy( dataDir )
% Liyun Tu Sep 25, 2014
%% Compute the geometry entropy using matrix as input.
% GEO_input.txt: a file which stored a matrix holding the p r u values of
% all the sreps. Each row is a feature, each column is a srep.
% dataDir = '/work/ltu/WorkSpace/Sep_3_newFeature/downSide/temp_data/';
% optionPNS = 1; small circle.


%%  Load data from .txt file, which generate by C++ containing matrix with  %%    
    GEO_data = load(strcat(dataDir,'/GEO_input.txt'));   

    if( isempty(GEO_data) )
        disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
    return ;
    else        
        disp(['S-rep Geometric features have been read from the folder: ' dataDir]) ;
    end     
   
    [rowNum,nSamples] = size(GEO_data);
    
    nSpokes = rowNum/7;    % rowNum=nSpokes*3(position) + nSpokes(radius) + nSpokes*3(direction)
    
    nPositions = 3 * nSpokes;
 
%% CPNS: Step 1 : Deal with hub Positions (PDM) %%    

    % position(i, j, k) =  j-th co-ordinate of i-th atoms position of the kth sample
    position = zeros( nSpokes, 3, nSamples ) ;
    
    for i = 1 : nSpokes,
        for j = 1 : 3,
            position( i, j, : ) = GEO_data( 3*(i-1)+j, : );
        end
    end     
    
    meanOfEachPDM = mean( position, 1 ) ;% 1X3XN matrix storing the COM of each sample    
   
    % translated positions centered at the mean of each sample    
    cposition = position - repmat( meanOfEachPDM, [nSpokes, 1, 1] ) ; 

    sscarryPDM = sqrt(sum(sum(cposition.^2))); %sum(sum(X)) return the sum of all the elements in X.
    
    % sphmatPDM = position part of the CPNS data (translated + scaled)
    % this is the input to the PNS program for the SHAPE part    
    sphmatPDM = zeros( nPositions, nSamples );
    spharrayPDM = zeros( nSpokes, 3, nSamples );

    for i = 1 : nSamples,
        spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
        sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nSpokes * 3, 1 );
    end

    % ------------------ Shape analysis via PNS ------------------  
     
    % compute PNS
    [ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100 );  %add default parameters values.
     
    disp('PNS computed for Shape PDMs') ;
    
    % Size of PDM's     
    sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
    sizePDM(:) = sscarryPDM(1, 1, :) ;    

    meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

    normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.
 
%% CPNS: Step 2 : Deal with atom radii (log radii are Euclidean variables) %%    

    logR = zeros( nSpokes, nSamples ) ;     % logR(i, j): i-th radius of j-th sample
    
    for i = 1 : nSpokes,        
        logR(i,:) = log(GEO_data(i+nPositions,:));
    end

    meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
    
    % r_bar(i, :): mean of i-th radius for all samples
    % used as a scaling factor for radii and spoke directions
    rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
    
    uScaleFactors =  zeros( 2 * nSpokes, nSamples ) ;
    for i = 1 : nSpokes,
        uScaleFactors( 2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
        uScaleFactors( 2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
    end    

    % shifted log r-i's
    RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  

%% CPNS: Step 3 : Deal with spoke directions (direction analysis) %% %%

        ZSpoke = zeros( 2 * nSpokes,  nSamples ) ;

        for r = 1 : nSpokes,

            a = GEO_data( nPositions+nSpokes+(r-1)*3+1 : nPositions+nSpokes+(r-1)*3+3,:);

            % renomalized the spoke direction to make x,y,z a unit vector.
            for count = 1 : nSamples,
                normalizeScale = sqrt(power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2));                
                a(1,count) = a(1,count) / normalizeScale;
                a(2,count) = a(2,count) / normalizeScale;
                a(3,count) = a(3,count) / normalizeScale;
            end  

            [Zr PNSr] = PNSmain(a,1); 
           
            ZSpoke( 2*(r-1) + 1 : 2*(r-1) + 2, : ) = Zr ;             
        end
        %disp(ZSpoke);    
      

%% CPNS: Step 4 : Construct composite Z matrix and perform PCA %% %%

    ZComp = [   meanSizePDM * ZShape ; ... 
                meanSizePDM * normalizedSizePDM ; ... 
                rScaleFactors .* RStar ; ... 
                uScaleFactors .* ZSpoke  ];
 
[ eigenVectors, pScores, eigenValues ]= princomp( ZComp', 'econ' ) ;


%% Calculate Geometry entorpy
    threshold = eigenValues/sum(eigenValues); 

    %eigenValues(threshold<1.0e-3) = [];
    eigenValues(threshold<0.01) = [];
    
    [k, cc] = size(eigenValues);
%     disp(rr);
    
    sumlognumda = sum(log(eigenValues));    
 
    geometryEntropy = k * 1.4189  + sumlognumda * 0.5;


       % each dlmwrite call opening and closing the file.
       dlmwrite(strcat(dataDir,'/geoEntropyResult.txt'), geometryEntropy, 'delimiter', '', 'precision',7);

end
