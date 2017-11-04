%function [ totalEntropy ] = upSpokeEntropy()
geoEntropyDir = '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/updateUV/';
regEntropyDir = '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/top/';
side =0;

addpath('/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/Correspondence/');

% Add these path to make sure the calculateCPNS_eigenvalue can works.
addpath('/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/functions/', '/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/ios/');
addpath('/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/', '/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/CPNS/');
addpath('/NIRAL/work/ltu/WorkSpace/pablo_matlab/shapeStat/CPNS/Sungkyu_additional_functions/subfunctions');
addpath('/NIRAL/tools/Matlab2011a/toolbox/stats/stats/');
    
    %% Regularity entropy
    % Calculate the srep area and volume entropy:
%     areavolumeDataDir = strcat(regEntropyDir, 'srepareavolume/');
%     areavolume = regularityEntropy(areavolumeDataDir,3);
%     
%     % Calculate the srep horizonal edges entropy:
%     horizonalDataDir = strcat(regEntropyDir, 'srephorizonal/');
%     horizonal = regularityEntropy(horizonalDataDir,2);
%     
%     % Calculate the srep vertical edges entropy:
%     verticalDataDir = strcat(regEntropyDir, 'srepvertical/');
%     vertical = regularityEntropy(verticalDataDir,2);
%     
%     regEntropy = areavolume + horizonal + vertical;
%     disp('The regularity entropy is: ');
%     disp(regEntropy);

    
    calculateCPNS( geoEntropyDir, 'test', 1 );
    disp('Finished CPNS');
    
    %% Geometry entropy, use small circle.
%     geoEntropy = calculateGeoEntropy(geoEntropyDir, 1, side);
%     disp('The geometry entropy is: ');
%     disp(geoEntropy);
    
    %% Total entropy
%     totalEntropy = geoEntropy - regEntropy;
%     disp('The total up entropy is: ');
%     disp(totalEntropy);
    
    %% Save the total entropy to a .txt file.
%     fid = fopen('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/updateUV/result.txt','w');
%     fprintf(fid,'%.25f\n',totalEntropy);
%     fclose(fid);