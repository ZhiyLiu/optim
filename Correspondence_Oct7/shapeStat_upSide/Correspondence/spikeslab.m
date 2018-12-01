%% This script use "Spike and slab variable selection" to classify objects
%  Author: Zhiyuan Liu
%  Date: 2018.11.2

close all; clear;
%% set up path
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

%% load in data
%% 1. collect all links from 10 cases, each case relates to 1 column
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

dimPerLink = 12;
num_rows = 12 * 68 * 3; % == 2448
num_cols = size(s_reps, 1);
feature_mat = [];

disp(['Loading features from all linking structures...']);

% load original group
fileName = 'middle_linking_structure.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = s_reps(i, :);
    feature_vector = loadLinkingStructure(s_rep_path, fileName);
    feature_mat = [feature_mat, feature_vector];

end

%load small shift group
fileName = 'diseased_linking_structure.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = s_reps(i, :);
    feature_vector = loadLinkingStructure(s_rep_path, fileName);
    feature_mat = [feature_mat, feature_vector];

end

% group features
feature_mat = reorderFeatures(feature_mat, 12);

% remove link to 
feature_mat(545:612, :) = [];

[d, n] = size(feature_mat);
y = ones(1, n);
y(1, 1:n/2) = -1;
output_data = [y; feature_mat];

output_file_name = '/playpen/workspace/demo_multiobject/spikeslab/multishape/in.csv';
csvwrite(output_file_name, output_data);
%% call spike and slab algorithm to output model
script_path = '/playpen/workspace/demo_multiobject/spikeslab/multishape/';
script_name = 'multishape.R';

script_full_path = fullfile(script_path, script_name);
% save data into file in the script path

s = system(script_full_path);
%!R CMD BATCH '/playpen/workspace/demo_multiobject/spikeslab/multishape/multishape.R'
if s ~= 0
    disp(['error when calling r script']);
end