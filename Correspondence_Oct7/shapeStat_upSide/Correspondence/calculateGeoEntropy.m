% side: 0, up spokes; 1, down spokes; 2, crest spokes.
%function [ geoEntropy] = calculateGeoEntropy( dataDir, outputFilename, optionPNS, side )
function [ geoEntropy] = calculateGeoEntropy( dataDir, optionPNS, side )
% dataDir = '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/moved_sreps/';
% outputFilename='output-plot_LV';
% optionPNS = 1;
% side =0;

if( nargin < 3 )
    optionPNS = 0 ;
end

    
%%  Get data from the data folder containing fitted s-rep models %%    
           
	[ CPNS_data CPNS_info ] = getCPNSDataMatrix( dataDir, side ) ;       

    if( isempty(CPNS_data) )
        disp('Error in reading CPNS data matrix! Program terminated without completion.') ;
    return ;
    else        
        disp(['CPNS Data has been read from the folder: ' dataDir]) ;
    end      
   
    %write into a file.
    %dlmwrite('/NIRAL/work/ltu/WorkSpace/CPNS_data.txt',CPNS_data,' ');

   
    nSamples = CPNS_info.nSamps ; 
    
    nRows = CPNS_info.nRows ;
    nCols = CPNS_info.nCols ;
    
    nTotalAtoms = nRows * nCols ;
    nEndAtoms = 2 * nRows + 2 * (nCols - 2) ; % same as 2 * nCols + 2 * (nRows - 2)
    nStdAtoms = nTotalAtoms - nEndAtoms ;        

    if(side==2)
        spokeNums =  nEndAtoms; 
    else
        spokeNums =  nEndAtoms + nStdAtoms; 
    end    
    
    nTotalPositions = 3 * spokeNums;
%     disp('------------------------the spoke number is: ');
%     disp(spokeNums);
 
%% CPNS: Step 1 : Deal with hub Positions (PDM) %%    

    % position(i, j, k) =  j-th co-ordinate of i-th atoms position of the kth sample
    position = zeros( spokeNums, 3, nSamples ) ;
    
    for i = 1 : spokeNums,
        for j = 1 : 3,                
            %position( i, j, : ) = CPNS_data( 3*(i-1)+j, : ) + noise(3*(i-1)+j, : );
            position( i, j, : ) = CPNS_data( 3*(i-1)+j, : );
        end
    end     
    
    meanOfEachPDM = mean( position, 1 ) ;% 1X3XN matrix storing the mean of PDM of each sample    
   
    % translated positions centered at the mean of each sample    
    cposition = position - repmat( meanOfEachPDM, [spokeNums, 1, 1] ) ;
    
    % cposition = position - repmat( meanOfCombinedPDM, [spokeNums, 1, nSamples] ) ; 

    % dibyendu
    % sscarryPDM = ( 1 X 1 X nSamples ) array for holding size of each PDM
    % sscarryPDM(k) = size of the PDM for k-th sample     
    % sscarryPDM is ( 1 X 1 X nSamples ) array, holding size of each srep
    sscarryPDM = sqrt(sum(sum(cposition.^2))); %sum(sum(X)) return the sum of all the elements in X.

    
    % sphmatPDM = position part of the CPNS data (translated + scaled)
    % this is the input to the PNS program for the SHAPE part    
    sphmatPDM = zeros( nTotalPositions, nSamples );
    spharrayPDM = zeros( spokeNums, 3, nSamples );

    for i = 1 : nSamples,
        spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
        sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', spokeNums * 3, 1 );
    end

    % ------------------ Shape analysis via PNS ------------------  
     
    % compute PNS
    [ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100 );  %add default parameters values.  
    %[ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100, true);
 
     
    disp('PNS computed for Shape PDMs') ;

    
    % Size of PDM's     
    sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
    sizePDM(:) = sscarryPDM(1, 1, :) ;    

    meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

    normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.
 
%% CPNS: Step 2 : Deal with atom radii (log radii are Euclidean variables) %%    

    logR = zeros( spokeNums, nSamples ) ;     % logR(i, j): i-th radius of j-th sample
    
    for i = 1 : spokeNums,        
        logR(i,:) = log(CPNS_data(i+nTotalPositions,:));
    end

    meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
    
    % r_bar(i, :): mean of i-th radius for all samples
    % used as a scaling factor for radii and spoke directions
    rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
    
    uScaleFactors =  zeros( 2 * spokeNums, nSamples ) ;
    for i = 1 : spokeNums,
        uScaleFactors( 2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
        uScaleFactors( 2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
    end
    
%             disp('The rScaleFactors is: ');
%     disp(rScaleFactors);
%     disp('The uScaleFactors is: ');
%     disp(uScaleFactors);
    % shifted log r-i's
    RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  

%% CPNS: Step 3 : Deal with spoke directions (direction analysis) %% %%

        ZSpoke = zeros( 2 * spokeNums,  nSamples ) ; 

        tOpStart = tic ;

        for r = 1 : spokeNums,

            a = CPNS_data( nTotalPositions+spokeNums+(r-1)*3+1 : nTotalPositions+spokeNums+(r-1)*3+3,:);
            %% ----------------------------------Begin, add noise to spoke direction. only used for test set.----------------
%             disp('The current Up spoke direction is: ');
%             disp(a);
            
            %% Re normalized the spoke direction to add noise to test set, to avoid no variablity in the population.            
            theta = 1.0e-5;
            % Generator settings, 'default' will generate the same values
            % when matlab restart. Same as rand('state', 0).
%            rng('default');%fixed noise. Any time generate the same random values.            

            du=rand(3, nSamples);
            %du = randn(3, nSamples); % use Normally distributed pseudorandom numbers. But, small sample number, maybe mean is not 0??
            
            a = a + theta*du;

            % renomalized the spoke direction to make x,y,z a unit vector.
            for count = 1 : nSamples,
                normalizeScale = sqrt(power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2));
                %normalizeScale = power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2);
                a(1,count) = a(1,count) / normalizeScale;
                a(2,count) = a(2,count) / normalizeScale;
                a(3,count) = a(3,count) / normalizeScale;
            end  
            %% ----------------------------------End, the up re-normalized part only used for test set.----------------
%             disp('The renomalized Up spoke direction is: ');
%             disp(a);
%             disp(power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2));


            [Zr PNSr] = PNSmain(a,1);

            tElapsed = toc( tOpStart ) ;
            tRemaining = ( tElapsed / r ) * ( spokeNums - r ) ;


                disp( ['Computing PNS for spoke ' num2str(r) ' of ' num2str(spokeNums)] ) ;
                if( tRemaining > 60 )
                    disp( ['Estimated time remaining is ' num2str(tRemaining/60) ' mins'] ) ;
                else
                    disp( ['Estimated time remaining is ' num2str(tRemaining) ' secs'] ) ;
                end
                if( r == spokeNums )
                    disp(['Total time elasped for PNS calculation of ' num2str(spokeNums) ' spokes = ' num2str(tElapsed/60) ' mins']) ;
                end 
           
             ZSpoke( 2*(r-1) + 1 : 2*(r-1) + 2, : ) = Zr ;             
        end
        %disp(ZSpoke);    
   

%% CPNS: Step 4 : Construct composite Z matrix and perform PCA %% %%

    ZComp = [   meanSizePDM * ZShape ; ... 
                meanSizePDM * normalizedSizePDM ; ... 
                rScaleFactors .* RStar ; ... 
                uScaleFactors .* ZSpoke  ];
            [r1,c1]=size(ZShape);
            [r2,c2]=size(normalizedSizePDM);
            [r3,c3]=size(RStar);
            [r4,c4]=size(ZSpoke);
            [r5,c5]=size(ZComp);
            
%             disp('-------------------------the size of ZShape is: ');
%             disp(r1);
%             disp(c1);
% 
%             disp('-------------------------the size of normalizedSizePDM is: ');
%             disp(r2);
%             disp(c2);
%             
%             disp('-------------------------the size of RStar is: ');
%             disp(r3);
%             disp(c3);
%             
%             disp('-------------------------the size of ZSpoke is: ');
%             disp(r4);
%             disp(c4);
%             
%             disp('-------------------------the size of ZComp is: ');
%             disp(r5);
%             disp(c5);

%% princomp method use corvariance matrix to calculate the eigenvalue. This
% output the same result as the following standard PCA method. Remind: the 
% eigen values gotten from princomp only keep large eigenvalues, cut-offs small ones. 
%    [ eigenVectors, pScores, eigenValues ]= princomp( ZComp', 'econ' ) ;
     
%     %% Save the covariance matrix to a .txt file.
%      fid = fopen('/NIRAL/work/ltu/WorkSpace/princomp-eig-3.txt','w');
%      fprintf(fid,'%.25f\n',eigenValues);
%      fclose(fid);
    
%% Standard PCA
% ZComp is a matrix with row as feature(dimesion), column as a case(sample).
[m,n] = size(ZComp);

% mean for each dimension(one row is one dimension), its a m by 1 matrix.
C = mean(ZComp,2);

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

     
     
%% PCA using correlation matrix instead of corvariance matrix
% % ZComp is a matrix with row as feature(dimesion), column as a case(sample).
% [m,n] = size(ZComp);
% % disp('The scaled matrix ZComp is: ');
% % disp(ZComp);
% 
% % mean for each dimension(one row is one dimension), its a m by 1 matrix.
% C = mean(ZComp,2);
% 
% % subtract off the row mean from A
% M = repmat(C,1,n); %repmat(C,1,n) means tile 1*n patch of C.
% centereddata = ZComp - M;
% 
% % delete the rows if its elements all zero.
% zindex = all(centereddata==0,2);
% %disp('The row index is: ');
% %disp(zindex);
% centereddata(zindex,:) = [];
% % disp('The centereddata is: ');
% % disp(centereddata);
% 
% 
% % calculate the covariance matrix
% trancentereddata = transpose(centereddata);
% 
% covariance_1 = 1/(n-1) * centereddata * trancentereddata;
% % make the covariance matrix symmetric.
% covariance = (covariance_1 + transpose(covariance_1))/2;
% %  disp('The covariance matrix for entropy is: ');
% %  disp(covariance);
% 
% % Convert covariance matrix to correlation matrix.
% correlation_1 = corrcov(covariance);
% correlation = (correlation_1 + transpose(correlation_1))/2;
% 
% % find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V)
% [PC, V] = eig(correlation);
% 
% % extract diagonal of matrix (eigenvalues) as vector
% V = diag(V);
%     
% % descent order
% [junk, rindices] = sort(-1*V);
% eigenvalue = V(rindices);
% [eRows, eCols] = size(eigenvalue);
% indexi = 1;
% for ind = 1:eRows
%     if(eigenvalue(ind,1)>1e-6)  %only keep the eigenvalue which is big than 1e-6.
%         EigenvalueCut(indexi,1) = eigenvalue(ind,1);
%         indexi = indexi +1;
%     end
% end
% disp('Using correlation matrix, The eigen value cuted is: ');
% disp(EigenvalueCut);
% 
% % Calculate the entropy: 
% geoEntropy = 0;
% for i = 1:(indexi-1)    %only use the eigenvalues big than 1e-6.
%     geoEntropy = geoEntropy + log(EigenvalueCut(i,1));
% end  
%  disp('The geometry entropy is: ');
%  disp(geoEntropy);



%% Calculate Geometry entorpy
    threshold = eigenValues/sum(eigenValues); 

    eigenValues(threshold<1.0e-3) = [];
    
%     % save the eigenvalues into file.
%     fid = fopen('/NIRAL/work/ltu/WorkSpace/eig-cutoff-4.txt','w');
%     fprintf(fid,'%.25f\n',eigenValues);
%     fclose(fid);

%     if (isempty(eigenValues))
%         disp('No eigenValue bigger than 1.0e-10, the Geometry entropy is set to -1000.');
%         geoEntropy = -1000;  % attention: is this value resonable?  
%     else
%         geoEntropy = sum(log(eigenValues));
%     end
 
    geoEntropy = sum(log(eigenValues));

