% Calculate the up or down spoke total entropy ( totalEntropy = GeometryEntropy - RegularityEntropy)
% We do not need to deside if its up or down spoke. We just put the correct
% m3d files in proper folder, the matlab code used to calculate regularity entropy do
% not need to know if it is up or down spokes, because the features in txt files used
% is the final result get from up or down spokes, there are the same.
% That's to say: just put the correct .txt files(features) in three
% seperate folder is ok.

% Input:
%   geoEntropyDir: the folder contains the m3d files you used to compute geometry entropy (e.g.
%   '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/aftermovewithdeltauv/')
%
%   regEntropyDir: the folder contains the features you used to compute
%   regularity entropy (e.g.'/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/'),
%   under this folder there should have three individual folder: srepareavolume, srephorizonal, srepvertical which 
%   store the three kinds of features seperately.
%
%   side: 0 up spoke; 1 down spoke.

function [ totalEntropy ] = spokeEntropy( geoEntropyDir, regEntropyDir, side )
    %% Geometry entropy, use small circle.
    geoEntropy = calculateCPNS_eigenvalue(geoEntropyDir, 1, side);

    %% Regularity entropy
    % Calculate the srep area and volume entropy:
    areavolumeDataDir = strcat(regEntropyDir, 'srepareavolume/');
    areavolume = regularityEntropy(areavolumeDataDir,3);
    
    % Calculate the srep horizonal edges entropy:
    horizonalDataDir = strcat(regEntropyDir, 'srephorizonal/');
    horizonal = regularityEntropy(horizonalDataDir,2);
    
    % Calculate the srep vertical edges entropy:
    verticalDataDir = strcat(regEntropyDir, 'srepvertical/');
    vertical = regularityEntropy(verticalDataDir,2);
    
    regEntropy = areavolume + horizonal + vertical;
    
    %% Total entropy
    totalEntropy = geoEntropy - regEntropy;