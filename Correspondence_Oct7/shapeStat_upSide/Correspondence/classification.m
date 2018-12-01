%%% This function implements classification of multi-object shape
%%% Two classes with Bayesian inference. First approach is to model
%%% p(x_i|C) is a Gaussian distribution N(mu_i, v_i)
%% set up environment
addpath('/playpen/workspace/demo_multiobject/AJIVE/SRVF')
addpath('/playpen/workspace/demo_multiobject/AJIVE/Smoothing')
addpath('/playpen/workspace/demo_multiobject/AJIVE/General')
addpath('/playpen/workspace/demo_multiobject/AJIVE/BatchAdjust')
addpath('/playpen/workspace/demo_multiobject/AJIVE/AJIVECode')
baseDir = '/playpen/software/shapestat/shapeStat';

addpath(    baseDir, ...
            [baseDir filesep 'ios'], ...
            [baseDir filesep 'functions'],...
            [baseDir filesep 'PGA' filesep 'output-plot'],...
            [baseDir filesep 'PGA'], ...
            [baseDir filesep 'CPNS'], ...
            [baseDir filesep 'CPNS' filesep 'Sungkyu_additional_functions'],...
            [baseDir filesep 'CPNS' filesep 'Sungkyu_additional_functions' filesep 'subfunctions']...        
        );
addpath( genpath( fullfile( pwd, 'BatchAdjust' ) ) );
addpath( genpath( fullfile( pwd, 'General' ) ) );
addpath( genpath( fullfile( pwd, 'Smoothing' ) ) );

%% load train data (47) and test data (10)
loadTrainData;
% loadTestData;

%% Training-Maximize likelihood estimate Gaussian parameters for each feature
disp(['Training...']);
% [idList, links_interested] = selectLinks(featureMatrix, dimPerLink, 9);
% cpnsLinkingFeature(links_interested, 12);

%% Test- Compute the likelihoods for given class 1 and class 2. Then select the larger one
% h1 = G1(x_1)*G1(x_2)...G1(x_p);
% h2 = G2(x_1)*G2(x_2)...G2(x_p);
% return max(h1,h2);


