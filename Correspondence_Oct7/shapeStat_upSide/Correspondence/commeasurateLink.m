function links = commeasurateLink(feature_mat, dimPerLink)

nSamples = size(feature_mat, 2);
nLinks = size(feature_mat, 1) / dimPerLink;

ZVector = zeros(nLinks, nSamples);

% deduce dimension for each link (one data point in high dimension) via PNS
pos_u = 3;
pt_mat = zeros(nLinks, 2, nSamples);
for i = 1 : nLinks
    start_pos_u = (i-1)*dimPerLink + pos_u;
    for j = 0 : 1
        pt_mat(i, j+1, :) = feature_mat(start_pos_u + j, :);
    end
end

meanOfXY = mean(pt_mat, 1);
cposition = pt_mat - repmat(meanOfXY, [nLinks, 1, 1]);
sscarryPDM = sqrt(sum(sum(cposition.^2)));

nPositions = 2 * nLinks;
sphmatPDM = zeros( nPositions, nSamples );
spharrayPDM = zeros( nLinks, 2, nSamples );

for i = 1 : nSamples
    spharrayPDM(:,:,i) = cposition(:,:,i)/sscarryPDM(1,1,i);
    sphmatPDM(:,i) = reshape(spharrayPDM(:,:,i)', nLinks * 2, 1 );
end

[ZVector PNSVector] = PNSmain( sphmatPDM);

% compute scale factor for link vector
sizePDM = zeros( 1, nSamples ) ;        % gamma(i)    
sizePDM(:) = sscarryPDM(1, 1, :) ; 

% map positive to (-inf, +inf)
meanSizePDM = exp(mean(log(sizePDM)));  % gamma_bar = GM of gamma(i)    

normalizedSizePDM = log( sizePDM / meanSizePDM ); % gamma(i)^* = log(gamma(i)/gamma_bar): normalized size variables. Euclidean variables scaling.

% deal with link length (always positive)
logR = zeros(nLinks, nSamples);
for i = 1: nLinks
    logR(i, :) = log(feature_mat((i-1)*dimPerLink + 1, :));
end

meanLogR = mean(logR, 2);

scaleFactor = repmat( exp( meanLogR ), [ 1, nSamples ] ) ;  
RStar = logR - repmat( meanLogR, 1, nSamples );  

links = [ meanSizePDM * ZVector; ...
          meanSizePDM * normalizedSizePDM; ...
          scaleFactor .* RStar...
    ];
end