%%%% This function is main entry of JIVE analysis of linking
%%%% structure
%%%% Author: Zhiyuan Liu
%%%% Date: Aug 6, 2018


%%% the structure of link is defined as purlin, namely: 
%%% (px, py, pz, ux, uy, uz, r, l, i, nx, ny, nz)

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

loadLinkData;

%% 2. select those links who always link to the same object
disp(['Selecting links to the same object among the population...']);
% idList is 0-based id
% links_interested are links who always link to same object among
% population
[idList, links_interested] = selectLinks(linkFeature, 5, 2);
saveIDList('/playpen/workspace/demo_multiobject/selected_links.txt', idList);

nSamples = size(links_interested,2);
% select spokes, dimension: 7
spokes_interested = [];
% filter out link_to feature from links_interested
% dimension of links reduce from 5 to 4
links_selected = [];
dimSpoke = 7; 
dimLink = 5;

remove_middle_to_left = 0; % remove thos links from middle to left
remove_middle_to_right = 2; % only study links who connect to left, remove link to right
remove_nothing = 100; % keep all links 
for i = 1:length(idList)
    % index for spokes
    index = idList(i) + 1; % for c++, I subtracted idList. Now need to add 1
    spokes_interested = [spokes_interested; shapeFeature((index-1) * dimSpoke + 1: index * dimSpoke, :)];
    links_selected = [links_selected; links_interested((i-1)*dimLink+1,:); links_interested((i-1)*dimLink+3:i*dimLink, :)];
%     link_to = links_interested((i-1)*dimLink+2, 1);
%     % the interested links already selected out, just filter out link_to
%     if link_to == remove_nothing
%         link_to_left = zeros(4, nSamples);
%         links_selected = [links_selected; link_to_left];
%     else
%         links_selected = [links_selected; links_interested((i-1)*dimLink+1,:); links_interested((i-1)*dimLink+3:i*dimLink, :)];
%     end
    
end
dimLink = 4;

spokes_selected = spokes_interested;

% Subtract spoke length from link length
index_spoke_length = 7;
for i = 1:length(idList)
    spoke_length = spokes_selected((i-1)*dimSpoke + index_spoke_length, :);
    links_selected((i-1)*dimLink + 1, :) = links_selected((i-1)*dimLink + 1, :) - spoke_length;
end

num_spokes = size(spokes_selected, 1) / dimSpoke;
%% 3. commensurate features
% spokes_normalized = commeasurateShape(spokes_interested, 7);
% links_normalized = commeasurateLink(links_interested, 5);

% Centering the data
mean_spokes_selected = repmat(mean(spokes_selected, 2), [1 nSamples]);
spokes_selected_centered = spokes_selected - mean_spokes_selected;
links_selected_centered = links_selected - repmat(mean(links_selected,2), [1 nSamples]);

% group features
spokes_reordered = reorderFeatures(spokes_selected_centered, dimSpoke);
links_reordered = reorderFeatures(links_selected_centered, dimLink);

% save to csv for spike and slab
saveFeatures2CSV(spokes_reordered, links_reordered);

%% 4. AJIVE analysis


close all;
data1 = spokes_reordered; %spokes_normalized;%spokes_interested;
data2 = links_reordered; %links_normalized; %links_interested;
AJIVEPreVisualMJ({data1},{[2 4 6]},113,{'Inner-object: [2 4 6]'}) ;
AJIVEPreVisualMJ({data2},{[2 4 6]},113,{'Inter-object: [2 4 6]'}) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%%%%%%%%%%%%%%%%%%%% AJIVE main for 2 blocks%%%%%%%%%%%%%%%%
caBlocks=[];
caBlocks{1} = data1 ;
caBlocks{2} = data2 ;
vRank = [4 6];
caDataSetNames = {'Inner-object' 'Inter-object'} ;
paramstruct = struct('dataname',{caDataSetNames}, ...
                           'nresample',1000, ...
                           'ioutput',[1 1 1 1 1 1 1 1 1]) ;
outstruct = AJIVEMainMJ(caBlocks,vRank,paramstruct);

%%%%%%%%%%%%%%%%%%%%%%% AJIVE for finer blocks %%%%%%%%%%%%%%%%%%%
% spoke: {p[3], u[3], r[1]}
% link:{l[1], n[3]}
% caBlocks=[];
% caBlocks{1} = data1(1:3, :) ; %skeletal points
% caBlocks{2} = data1(4:6, :) ; % spoke dir
% caBlocks{3} = data1(7, :);    % spoke length
% caBlocks{4} = data2(1, :);    % link length
% caBlocks{5} = data2(2:end, :);% link vector
% vRank = [3 2 1 1 2];
% caDataSetNames = {'skeletal points' 'spoke dir' 'spoke length' 'link length' 'link vector'} ;
% paramstruct = struct('dataname',{caDataSetNames}, ...
%                            'nresample',1000, ...
%                            'ioutput',[1 1 1 1 1 1 1 1 1]) ;
% outstruct = AJIVEMainMJ(caBlocks,vRank,paramstruct);
% 
%% 5. draw loadings
figure; hold on;
subplot(3,2,1);
cns_spokes = outstruct.CNSloading{1,1};
feature_name = {'px','py','pz','ux','uy','uz','r'};
%barcolor=[0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
barcolor = ['b', 'm', 'c','r', 'g', 'y','k'];
titlestr = {'CNS Loadings of block1'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% CNS loading of block 2
subplot(3,2,2);
cns_spokes = outstruct.CNSloading{1,2};
feature_name = {'l','nx','ny','nz'};
barcolor = ['b', 'm', 'c','r'];
titlestr = {'CNS Loadings of block2'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% BSS joint loading of block 1
subplot(3, 2, 3);
cns_spokes = outstruct.BSSjointLoading{1,1};
feature_name = {'px','py','pz','ux','uy','uz','r'};
%barcolor=[0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
barcolor = ['b', 'm', 'c','r', 'g', 'y','k'];
titlestr = {'BSS joint Loadings of block1'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% BSS joint loading of block 2
subplot(3, 2, 4);
cns_spokes = outstruct.BSSjointLoading{1,2};
feature_name = {'l','nx','ny','nz'};
barcolor = ['b', 'm', 'c','r'];
titlestr = {'BSS joint Loadings of block2'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% BSS indiv loading of block 1
subplot(3, 2, 5);
cns_spokes = outstruct.BSSindivLoading{1,1};
feature_name = {'px','py','pz','ux','uy','uz','r'};
%barcolor=[0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
barcolor = ['b', 'm', 'c','r', 'g', 'y','k'];
titlestr = {'BSS individual Loadings of block1'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% BSS joint loading of block 2
subplot(3, 2, 6);
cns_spokes = outstruct.BSSindivLoading{1,2};
feature_name = {'l','nx','ny','nz'};
barcolor = ['b', 'm', 'c','r'];
titlestr = {'BSS individual Loadings of block2'};
drawBarGraph(cns_spokes, num_spokes, feature_name, barcolor, titlestr);

% draw row curves and cols curves
figure; hold on;
jc = outstruct.MatrixJoint{1,1};
[d n] = size(jc);
clrm = RainbowColorsQY(n);
colormap(clrm);
% cols of JC

for i = 1:nSamples
    plot(jc(:, i));
end
title('Columns as curves in joint component of inter-obj');
%draw JC scores
% each point has one color
% mcolor = [0 0 1;...
%     0 0.5 1;...
%     0 1 0;...
%     0.5 1 0.5;...
%     0 1 1;...
%     0.2 0.6 1;...
%     1 0 0;...
%     1 0 1;...
%     1 1 0;...
%     1 1 0.3;...
%     1 1 1];
close all;
% JC 1 scores
figure; hold on;
draw_data = outstruct.MatrixJoint{1,1};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'Joint matrix for data block 1'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);

% JC 2 scores
figure; hold on;
draw_data = outstruct.MatrixJoint{1,2};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'Joint matrix for data block 2'}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);

% JC 1 loadings
figure; hold on;
draw_data = outstruct.MatrixJoint{1,1};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'Joint matrix for data block 1'}}
drawScatPlot(draw_data', mcolor, legendstr, markerstr, titlestr);
% JC 2 loadings
figure; hold on;
draw_data = outstruct.MatrixJoint{1,2};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'Joint matrix for data block 2'}}
drawScatPlot(draw_data', mcolor, legendstr, markerstr, titlestr);

%indiv scores
figure; hold on;
draw_data = outstruct.MatrixIndiv{1,1};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{''}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);
%indiv scores
figure; hold on;
draw_data = outstruct.MatrixIndiv{1,2};
legendstr = {{'orginal 10 cases'}};

mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{''}}
drawScatPlot(draw_data, mcolor, legendstr, markerstr, titlestr);
