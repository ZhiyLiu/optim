%% 1. collect all links from 10 cases, each case relates to 1 column
data_path = '/playpen/workspace/proj_multiobject/';
test_case = ['889945'; '895096'; '910593'; '915717';'931663';'938659';'950194';'961753';'967950';'974849'];

%% 2. load class 1
num_cols = size(train_case, 1);
fileName = 'middle_linking_structure.m3d';
featureMatrix = [];
disp(['Loading features from all linking structures...']);
for i = 1 : size(train_case, 1)
    s_rep_path = fullfile(data_path, train_case(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    featureVector = [shape; link];
    featureMatrix = [featureMatrix, featureVector];
end

%% 3. load class 2
fileName = 'diseased_linking_2.m3d';
for i = 1 : size(train_case, 1)
    s_rep_path = fullfile(data_path, train_case(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    featureVector = [shape; link];
    featureMatrix = [featureMatrix, featureVector];
end