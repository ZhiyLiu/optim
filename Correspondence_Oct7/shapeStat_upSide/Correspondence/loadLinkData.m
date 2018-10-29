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

legend_str = {{'180215','223443','225196','402809','549967','622437','721965','730231','802785','867237',...
    '180215_diseased','223443_diseased','225196_diseased','402809_diseased',...
    '549967_diseased','622437_diseased','721965_diseased','730231_diseased',...
    '802785_diseased','867237_diseased'}};
% 12-dimensional(purlz+linkTo) links and 68 links per object
% each column consists of feature vector per configuration(3 objects)

num_rows = 12 * 68 * 3; % == 2448
num_cols = size(s_reps, 1);
fileName = 'middle_linking_structure.m3d';
shapeFeature = [];
linkFeature = [];
disp(['Loading features from all linking structures...']);
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(s_reps(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    shapeFeature = [shapeFeature, shape];
    linkFeature = [linkFeature, link];
end

% load 1st disease group
fileName = 'diseased_linking_structure.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(s_reps(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    shapeFeature = [shapeFeature, shape];
    linkFeature = [linkFeature, link];
end


% load 2nd 'diseased' linkinking structure
% fileName = 'diseased_linking_2.m3d';
% for i = 1 : size(s_reps, 1)
%     s_rep_path = fullfile(s_reps(i, :),fileName);
%     [shape link] = loadJIVE(s_rep_path);
%     shapeFeature = [shapeFeature, shape];
%     linkFeature = [linkFeature, link];
% end

%% Study links between 2 objects
