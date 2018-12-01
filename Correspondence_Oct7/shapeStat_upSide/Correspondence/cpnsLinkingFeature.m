%%% This function normalize feature vectors of linking structure
%%% each column of feature_mat contains links from 3 object in this case
%%% dimPerLink is predefined value: the dimension of one single link
%%% the structure of link is defined as purlin, namely: 
%%% input (px, py, pz, ux, uy, uz, r, l, i, nx, ny, nz)
function cpnsLinkingFeature(feature_mat, dimPerLink)

%% Normalize position
nSamples = size(feature_mat, 2);
nLinks = size(feature_mat, 1) / dimPerLink;

% in 2D, all z value should be the same, no variation, dismiss
position_mat = ones(nLinks, 2, nSamples);
for i = 1: nLinks
    for j = 1:2
        position_mat(i, j, :) = feature_mat((i-1)*dimPerLink + j, :);
    end
end

% centering
meanOfXY = mean(position_mat, 1); % column 1 = x coordinate; 2=y; 3=z
meanOfCombinedPDM = mean( mean( position_mat, 1 ), 3 ) ;
cposition = position_mat - repmat(meanOfXY, [nLinks, 1, 1]); %centering

sscarryPDM = sqrt(sum(sum(cposition.^2))); % for each sample, compute the standard deviation


% divide by variance
% sphmatPDM = position part of the CPNS data (translated + scaled)
% this is the input to the PNS program for the SHAPE part    
nPositions = 2 * nLinks;
sphmatPDM = zeros( nPositions, nSamples );
spharrayPDM = zeros( nLinks, 2, nSamples );

for i = 1 : nSamples,
    spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
    sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nPositions, 1 );
end

% ------------------ Shape analysis via PNS ------------------  

% compute PNS
[ZShape PNSShape] = PNSmain( sphmatPDM, 0, 0.05, 100 );  %add default parameters values.

disp('PNS computed for skeletal position of links') ;
%% scale factor ( deal with standard deviation)
sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
sizePDM(:) = sscarryPDM(1, 1, :) ;    

% map positive to (-inf, +inf)
meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.

%% Euclideanize spoke direction
PNSSpoke = cell( nLinks, 1 ) ;

%ZSpoke = zeros(2 * nLinks,  nSamples ) ; %3D
ZSpoke = zeros(nLinks, nSamples); %2D
pos_u = 4;

for r = 1 : nLinks,

    start_pos_u = (r-1)*dimPerLink + pos_u;
    a = feature_mat( start_pos_u : start_pos_u + 1, :);

    % renomalized the spoke direction to make x,y,z a unit vector.
%     for count = 1 : nSamples,
%         normalizeScale = sqrt(power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2));                
%         a(1,count) = a(1,count) / normalizeScale;
%         a(2,count) = a(2,count) / normalizeScale;
%         a(3,count) = a(3,count) / normalizeScale;
%     end  

    [Zr PNSr] = PNSmain(a,0, 0.05, 100, true); 

    PNSSpoke{r} = PNSr;

    %ZSpoke( 2*(r-1) + 1 : 2*(r-1) + 2, : ) = Zr ;  
    ZSpoke( r,: ) = Zr(1,:) ;  
end

%% Euclideanize spoke length
logR = zeros( nLinks, nSamples ) ;     % logR(i, j): i-th radius of j-th sample

pos_r = 7; % the position of R in link structure
for i = 1 : nLinks,        
    logR(i,:) = log(feature_mat((i-1)* dimPerLink+pos_r,:));
end

meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples
meanRs = exp(meanLogR);

% r_bar(i, :): mean of i-th radius for all samples
% used as a scaling factor for radii and spoke directions
rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  

% uScaleFactors =  zeros(2 * nLinks, nSamples ) ;
% for i = 1 : nLinks,
%     uScaleFactors(2*(i-1) + 1, : ) = rScaleFactors( i, : ) ;
%     uScaleFactors(2*(i-1) + 2, : ) = rScaleFactors( i, : ) ;
% end   

% shifted log r-i's
RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  

%% Euclideanize link length
logL = zeros(nLinks, nSamples);
pos_l = 8;
for i = 1: nLinks
    logL(i,:) = log(feature_mat((i-1)*dimPerLink+pos_l, :));
end

meanLogL = mean(logL, 2);
LStar = logL - repmat( meanLogL, 1, nSamples);
lScaleFactors = repmat(exp(meanLogL), [1, nSamples]);

%% Normalize link point
link_position_mat = ones(nLinks, 2, nSamples);
pos_z = 10;
for i = 1: nLinks
    for j = 1:2
        
        link_position_mat(i, j, :) = feature_mat((i-1)*dimPerLink + j-1 + pos_z, :);
    end
end

% centering
meanOfXY = mean(link_position_mat, 1); % column 1 = x coordinate; 2=y; 3=z
cposition = link_position_mat - repmat(meanOfXY, [nLinks, 1, 1]); %centering

sscarryPDM = sqrt(sum(sum(cposition.^2))); % for each sample, compute the standard deviation


% divide by variance
% sphmatPDM = position part of the CPNS data (translated + scaled)
% this is the input to the PNS program for the SHAPE part    
nPositions = 2 * nLinks;
sphmatPDM = zeros( nPositions, nSamples );
spharrayPDM = zeros( nLinks, 2, nSamples );

for i = 1 : nSamples,
    spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
    sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nLinks * 2, 1 );
end

% ------------------ Shape analysis via PNS ------------------  

% compute PNS
[ZLinkShape PNSLinkShape] = PNSmain( sphmatPDM, 0, 0.05, 100 );  %add default parameters values.

disp('PNS computed for skeletal position of links') ;
%% scale factor ( deal with standard deviation)
sizeLink = zeros( 1, nSamples ) ;        % gamma(i)    
sizeLink(:) = sscarryPDM(1, 1, :) ;    

% map positive to (-inf, +inf)
meanSizeLink = exp(mean(log(sizeLink)));  % gamma_bar = GM of gamma(i)    

normalizedSizeLink = log( sizeLink/ meanSizeLink); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.


%% Now is ready to construct composite Z matrix and perform PCA
ZComp = [meanSizePDM * ZShape; ...
          meanSizePDM * normalizedSizePDM; ...
          rScaleFactors .* ZSpoke; ...
          rScaleFactors .* RStar;...
          lScaleFactors .* LStar; ...
          meanSizeLink * ZLinkShape; ...
          meanSizeLink * normalizedSizeLink; ...
          ];

%% 4. apply PCA and render

% The maximum no. of CPNS modes to be written to the CPNS mean file
MAX_CPNS_MODES = 15 ; 

[ eigenVectors, pScores, eigenValues ]= princomp( ZComp', 'econ' ) ;

nEigenModesToWrite = min( MAX_CPNS_MODES, length(eigenValues) ) ;
% 
% if( subtractNoisyVariance )        
%     if( nEigenModesToWrite < length( eigenValues ) )
%         eigenValToSubtract = mean( eigenValues(nEigenModesToWrite+1:end) ) ;
%         eigenValues = eigenValues - eigenValToSubtract ;
%         eigenValues( eigenValues < 0 ) = 0.0 ;
%     end        
% end


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

% [pathstr, name, ext] = fileparts(outputFilename) ;
% 
% imgEigenValsFilename = fullfile(pathstr, [name '_CPNS_eVals.jpg']) ;    
% saveas( hFig, imgEigenValsFilename ) ;
% disp(['Eigenvalues image written to: ' imgEigenValsFilename ]) ;
    
%% CPNS: Step 5: Transform mean back from E-comp (eucliden) space to S (SRep) space %% 
% 1. Convert from EComp to S-space

CPNSMeanScores = zeros( size(eigenVectors, 1), 1 ) ;     % = z vector    

CPNS_info = struct('nLinks', nLinks, 'nSamples', nSamples);

%linkData = convertECompToLinkData( CPNS_info, CPNSMeanScores, meanSizePDM, meanOfCombinedPDM, meanRs, PNSShape, PNSSpoke ) ;     


end