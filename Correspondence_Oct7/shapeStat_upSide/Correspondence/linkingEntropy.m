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
data_path = '/playpen/workspace/demo_multiobject/';
s_reps = ['180215';
    '223443';
    '225196';
    '402809';
    '549967';
    '622437';
    '721965';
    '730231';
    '802785';
    '867237'];

% s_reps = ['107524'; '141335'; '145680'; '180215';'223443';'225196';'236316';'252848';'315149';'323623';
%           '339961'; '359973'; '360222'; '365557';'402809';'429591';'482642';'490181';'513330';'549587';
%           '549967'; '589344'; '618672'; '622437';'631202';'638999';'641078';'660372';'691501';'699209';
%           '713824'; '721965'; '730231'; '755016';'792210';'793001';'802785';'810313';'812952';'813346';
%           '889945'; '895096'; '910593'; '915717';'931663';'938659';'950194';'961753';'967950';'974849';
%           '822794'; '833557'; '841812'; '842426';'867237';'876501';'879873'];
% 12-dimensional(purlz+linkTo) links and 68 links per object
% each column consists of feature vector per configuration(3 objects)
dimPerLink = 12;
num_rows = 12 * 68 * 3; % == 2448
num_cols = size(s_reps, 1);
feature_mat = [];

disp(['Loading features from all linking structures...']);
fileName = 'middle_linking_structure.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(data_path,s_reps(i, :));
    feature_vector = loadLinkingStructure(s_rep_path, fileName);
    feature_mat = [feature_mat, feature_vector];

end

% load 1st 'diseased' linkinking structure
% fileName = 'diseased_linking_structure.m3d';
% for i = 1 : size(s_reps, 1)
%     s_rep_path = fullfile(data_path,s_reps(i, :));
%     feature_vector = loadLinkingStructure(s_rep_path,fileName);
%     feature_mat = [feature_mat, feature_vector];
% end

% load 2nd 'diseased' linkinking structure
fileName = 'diseased_linking_2.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(data_path,s_reps(i, :));
    feature_vector = loadLinkingStructure(s_rep_path,fileName);
    feature_mat = [feature_mat, feature_vector];
end

%% 2. select those links who always link to the same object
disp(['Selecting links to the same object among the population...']);
% idList is 0-based id
[idList, links_interested] = selectLinks(feature_mat, dimPerLink, 9);
saveIDList('/playpen/workspace/demo_multiobject/selected_links.txt', idList);


% Subtract spoke length from link length
% index_spoke_length = 7;
% index_link_length = 8;
% for i = 1:length(idList)
%     links_interested((i-1)*dimPerLink + index_link_length, :) = links_interested((i-1)*dimPerLink + index_link_length, :) - links_interested((i-1)*dimPerLink + index_spoke_length, :);
% end

%% 3. euclideanize feature matrix
disp(['Commensurating and composing Z matrix...']);
cpnsLinkingFeature(links_interested, dimPerLink);

% show eigen modes of each block of the matrix by commenting out some code
%cpnsShapeFeature(links_interested, dimPerLink);

