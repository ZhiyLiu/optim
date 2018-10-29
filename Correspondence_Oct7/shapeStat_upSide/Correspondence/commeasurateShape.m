function shape = commeasurateShape(feature_mat, dimPerLink)
shape = [];
%% Normalize position
nSamples = size(feature_mat, 2);
nLinks = size(feature_mat, 1) / dimPerLink;

position_mat = ones(nLinks, 3, nSamples);
for i = 1: nLinks
    for j = 1:3
        position_mat(i, j, :) = feature_mat((i-1)*dimPerLink + j, :);
    end
end

% centering
meanOfXYZ = mean(position_mat, 1); % column 1 = x coordinate; 2=y; 3=z
cposition = position_mat - repmat(meanOfXYZ, [nLinks, 1, 1]); %centering

sscarryPDM = sqrt(sum(sum(cposition.^2))); % for each sample, compute the standard deviation


% divide by variance
% sphmatPDM = position part of the CPNS data (translated + scaled)
% this is the input to the PNS program for the SHAPE part    
nPositions = 3 * nLinks;
sphmatPDM = zeros( nPositions, nSamples );
spharrayPDM = zeros( nLinks, 3, nSamples );

for i = 1 : nSamples,
    spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
    sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nLinks * 3, 1 );
end

% ------------------ Shape analysis via PNS ------------------  

% compute PNS
[ZShape PNSShape] = PNSmain( sphmatPDM, 1, 0.05, 100 );  %add default parameters values.

disp('PNS computed for skeletal position of links') ;
%% scale factor ( deal with standard deviation)
sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
sizePDM(:) = sscarryPDM(1, 1, :) ;    

% map positive to (-inf, +inf)
meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.

%% Euclideanize spoke direction
ZSpoke = zeros(nLinks,  nSamples ) ;
pos_u = 4;

for r = 1 : nLinks,

    start_pos_u = (r-1)*dimPerLink + pos_u;
    a = feature_mat( start_pos_u : start_pos_u + 2, :);

    % renomalized the spoke direction to make x,y,z a unit vector.
%     for count = 1 : nSamples,
%         normalizeScale = sqrt(power(a(1,count),2) + power(a(2,count),2) +power(a(3,count),2));                
%         a(1,count) = a(1,count) / normalizeScale;
%         a(2,count) = a(2,count) / normalizeScale;
%         a(3,count) = a(3,count) / normalizeScale;
%     end  

    [Zr PNSr] = PNSmain(a,1); 

    ZSpoke( r, : ) = Zr(1,:) ;             
end

%% Euclideanize spoke length
logR = zeros( nLinks, nSamples ) ;     % logR(i, j): i-th radius of j-th sample

pos_r = 7; % the position of R in link structure
for i = 1 : nLinks,        
    logR(i,:) = log(feature_mat((i-1)* dimPerLink+pos_r,:));
end

meanLogR = mean(logR,2);                    % logR_bar(i): mean of i-th radius along all samples

% r_bar(i, :): mean of i-th radius for all samples
% used as a scaling factor for radii and spoke directions
rScaleFactors = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  

uScaleFactors =  zeros(nLinks, nSamples ) ;
for i = 1 : nLinks,
    uScaleFactors(i, : ) = rScaleFactors( i, : ) ;
    %uScaleFactors(2*(r-1) + 1, : ) = rScaleFactors( i, : ) ;
    %uScaleFactors(2*(r-1) + 2, : ) = rScaleFactors( i, : ) ;
end    

% shifted log r-i's
RStar = logR - repmat( meanLogR, 1, nSamples );       % R*(i, j) = logR(i, j) - logR_bar(i)  
shape = [meanSizePDM * ZShape; ...
          meanSizePDM * normalizedSizePDM; ...
          rScaleFactors .* RStar;...
          uScaleFactors .* ZSpoke; ...
          ];
end