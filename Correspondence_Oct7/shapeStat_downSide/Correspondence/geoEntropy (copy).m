%function [ geometryEntropy ] = geoEntropy( dataDir )
% Liyun Tu Sep 9, 2014

%% Compute the geometry entropy using matrix as input, not read from m3d file.
% the srep info was put into 3 files.
% GEO_hubposition.txt: 

% CPNS_input_matrix.txt: a file which stored a matrix holding the p r u values of all the sampled

% dataDir = '/work/ltu/WorkSpace/Sep_3_newFeature/downSide/temp_data/';
% optionPNS = 1; small circle.

% returns the full path and name of the file in which the call occurs, not including the filename extension.
% p = /work/ltu/WorkSpace/project/shapeStat/Correspondence/geoEntropy
p = mfilename('fullpath');

%[pathstr,name,ext] = fileparts(p) 

% %parent_path = /work/ltu/WorkSpace/project/shapeStat/Correspondence
%parent_path = pathstr(1:end-14); % /work/ltu/WorkSpace/project/shapeStat/
parent_path = p(1:end-25);

addpath(strcat(parent_path, 'functions/'));
addpath(strcat(parent_path, 'ios/'));
addpath(parent_path);
addpath(strcat(parent_path, 'CPNS/'));
addpath(strcat(parent_path, 'CPNS/Sungkyu_additional_functions/subfunctions/'));


%%  Load data from .txt file, which generate by C++ containing matrix with  %%  
    GEO_hubpos = load(strcat(dataDir,'/GEO_hubposition.txt'));

    if( isempty(GEO_hubpos) )
        disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
    return ;
    else        
        disp(['CPNS Data has been read from the folder: ' dataDir]) ;
    end      

    [nSamples, colNum] = size(GEO_hubpos);

    nSpokes = colNum/3;
 
%% CPNS: Step 1 : Deal with hub Positions (PDM) %%    

    % position(i, j, k) =  j-th co-ordinate of i-th atoms position of the kth sample
    position = zeros( nSpokes, 3, nSamples ) ;
    
    for i = 1 : nSpokes,
        for j = 1 : 3,
            position( i, j, : ) = GEO_hubpos(:, 3*(i-1)+j);
        end
    end     
    
    meanOfEachPDM = mean( position, 1 ) ;% 1X3XN matrix storing the mean of PDM of each sample    
   
    % translated positions centered at the mean of each sample    
    cposition = position - repmat( meanOfEachPDM, [nSpokes, 1, 1] ) ;

    sscarryPDM = sqrt(sum(sum(cposition.^2))); %sum(sum(X)) return the sum of all the elements in X.

    
    % sphmatPDM = position part of the CPNS data (translated + scaled)
    % this is the input to the PNS program for the SHAPE part    
    sphmatPDM = zeros( colNum, nSamples );
    spharrayPDM = zeros( nSpokes, 3, nSamples );

    for i = 1 : nSamples,
        spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
        sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nSpokes * 3, 1 );
    end


    % compute PNS
    [ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100 );  %add default parameters values.
     
    disp('PNS computed for Shape PDMs') ;
    
    % Size of PDM's     
    sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
    sizePDM(:) = sscarryPDM(1, 1, :) ;    

    meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

    normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.
 
%% CPNS: Step 2 : Deal with atom radii (log radii are Euclidean variables) %%  
    GEO_radius = load(strcat(dataDir,'/GEO_radius.txt'));

    if( isempty(GEO_radius) )
        disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
    return ;
    else        
        disp(['CPNS Data has been read from the folder: ' dataDir]) ;
    end

    logR = zeros( nSpokes, nSamples ) ;     % logR(i, j): i-th radius of j-th sample
    
    for i = 1 : nSpokes,        
        logR(i,:) = log(GEO_radius(i, :));
    end
    
    meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
    
    % r_bar(i, :): mean of i-th radius for all samples
    % scaling factor for radii and spoke directions
    rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
    
    uScaleFactors =  zeros( 2 * nSpokes, nSamples ) ;
    for i = 1 : nSpokes,
        uScaleFactors( 2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
        uScaleFactors( 2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
    end    

    % shifted log r-i's
    RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  

%% CPNS: Step 3 : Deal with spoke directions (direction analysis) %% %%
    GEO_dir = load(strcat(dataDir,'/GEO_direction.txt'));
    
    if( isempty(GEO_dir) )
        disp('Error in reading GEO_dir matrix! Program terminated without completion.') ;
    return ;
    else        
        disp(['CPNS Data has been read from the folder: ' dataDir]) ;
    end

        ZSpoke = zeros( 2 * nSpokes,  nSamples ) ; 

        for r = 1 : nSpokes,

            % for each spoke across the training set
            a = GEO_dir(:, (r-1)*3+1 : (r-1)*3+3);
            
            a = a';
            
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
  
disp('PNS computed for Direction') ; 
   

%% CPNS: Step 4 : Construct composite Z matrix and perform PCA %% %%

    ZComp = [   meanSizePDM * ZShape ; ... 
                meanSizePDM * normalizedSizePDM ; ... 
                rScaleFactors .* RStar ; ... 
                uScaleFactors .* ZSpoke  ];
 

% Delete the 0 rows in ZComp. (Because 3*39-4 > 30, the PDM sphere is 29-dimension, the rest of 3*39-4-29 rows is all 0)
zindex = all(ZComp < 0.000001, 2); % the rows indexs with all its elements zero.
ZComp(zindex,:) = [];

% save ZComp to a .txt file
%dlmwrite('/NIRAL/work/ltu/WorkSpace/CorrespondenEvaluation/beforeCorr/ZComp.txt',ZComp,' ');

    
%% Standard PCA
% ZComp is a matrix with row as feature(dimesion), column as a case(sample).
[m,n] = size(ZComp);

%              disp('-------------------------the size of ZComp is: ');
%              disp(ZComp);

% mean for each dimension(one row is one dimension), its a m by 1 matrix.
C = mean(ZComp,2);
disp('---c should be 0');
disp(C);

% subtract off the row mean from A
M = repmat(C,1,n); %repmat(C,1,n) means tile 1*n patch of C.
centereddata = ZComp - M;

% calculate the covariance matrix
trancentereddata = transpose(centereddata);
covariance_1 = 1/(n-1) * centereddata * trancentereddata;
% make the covariance matrix symmetric.
covariance = (covariance_1 + transpose(covariance_1))/2;

% Save the covariance matrix to a .txt file.
% dlmwrite('/NIRAL/work/ltu/WorkSpace/covMatrix-7.txt',covariance,' ');

[PC, V] = eig(covariance);

% extract diagonal of matrix (eigenvalues) as vector
V = diag(V);

% sort the variances in decreasing order
% rindices is the original index of element in V, use this we can get the
% element from V in decreasing order.
[junk, rindices] = sort(-1*V);
eigenValues = V(rindices);


%% Calculate Geometry entorpy
    threshold = eigenValues/sum(eigenValues); 

    eigenValues(threshold<1.0e-3) = [];
    
    [k, cc] = size(eigenValues);
%     disp(rr);
    
    sumlognumda = sum(log(eigenValues));    
 
    geometryEntropy = k * 1.4189  + sumlognumda * 0.5;

%end

    %% Save the total entropy to a .txt file.
     fid = fopen(strcat(dataDir,'/geoEntropyResult.txt'),'w');
     fprintf(fid,'%.25f\n',geometryEntropy);
     fclose(fid);

quit;
