%%%% This function is main entry of analyzing linking structure with CPNS
%%%% 
%%%% Author: Zhiyuan Liu
%%%% Date: Aug 6, 2018


%%% the structure of link is defined as purlin, namely: 
%%% (px, py, pz, ux, uy, uz, r, l, i, nx, ny, nz)
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

%% 1. collect all links from 10 cases, each case relates to 1 column
clear all;
close all;
s_reps = ['/playpen/workspace/demo_multiobject/180215';
    '/playpen/workspace/demo_multiobject/223443';
    '/playpen/workspace/demo_multiobject/225196';
    '/playpen/workspace/demo_multiobject/402809';
    '/playpen/workspace/demo_multiobject/549967';
    '/playpen/workspace/demo_multiobject/622437';
    '/playpen/workspace/demo_multiobject/721965';
    '/playpen/workspace/demo_multiobject/730231';
    '/playpen/workspace/demo_multiobject/802785';
    '/playpen/workspace/demo_multiobject/867237'];

% 12-dimensional(purlz+linkTo) links and 68 links per object
% each column consists of feature vector per configuration(3 objects)
dimPerLink = 12;
num_rows = 12 * 68 * 3; % == 2448
num_cols = size(s_reps, 1);
feature_mat = [];

disp(['Loading features from all linking structures...']);
fileName = 'middle_linking_structure.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = s_reps(i, :);
    feature_vector = loadLinkingStructure(s_rep_path, fileName);
    feature_mat = [feature_mat, feature_vector];

end

% load 1st 'diseased' linkinking structure
% fileName = 'diseased_linking_structure.m3d';
% for i = 1 : size(s_reps, 1)
%     s_rep_path = s_reps(i, :);
%     feature_vector = loadLinkingStructure(s_rep_path,fileName);
%     feature_mat = [feature_mat, feature_vector];
% end

% load 2nd 'diseased' linkinking structure
% fileName = 'diseased_linking_2.m3d';
% for i = 1 : size(s_reps, 1)
%     s_rep_path = s_reps(i, :);
%     feature_vector = loadLinkingStructure(s_rep_path,fileName);
%     feature_mat = [feature_mat, feature_vector];
% end

%% 2. select those links who always link to the same object
disp(['Selecting links to the same object among the population...']);
% idList is 0-based id
[idList, links_interested] = selectLinks(feature_mat, dimPerLink, 9);
saveIDList('/playpen/workspace/demo_multiobject/selected_links.txt', idList);


% Subtract spoke length from link length
index_spoke_length = 7;
index_link_length = 8;
for i = 1:length(idList)
    links_interested((i-1)*dimPerLink + index_link_length, :) = links_interested((i-1)*dimPerLink + index_link_length, :) - links_interested((i-1)*dimPerLink + index_spoke_length, :);
end

%% 3. euclideanize feature matrix
disp(['Commensurating and composing Z matrix...']);
cpnsLinkingFeature(links_interested, dimPerLink);

% show eigen modes of each block of the matrix by commenting out some code
%cpnsShapeFeature(links_interested, dimPerLink);

