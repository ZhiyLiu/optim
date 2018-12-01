%% 1. collect all links from 47 cases, each case relates to 1 column
data_path = '/playpen/workspace/proj_multiobject/';
train_case = ['107524'; '141335'; '145680'; '180215';'223443';'225196';'236316';'252848';'315149';'323623';
          '339961'; '359973'; '360222'; '365557';'402809';'429591';'482642';'490181';'513330';'549587';
          '549967'; '589344'; '618672'; '622437';'631202';'638999';'641078';'660372';'691501';'699209';
          '713824'; '721965'; '730231'; '755016';'792210';'793001';'802785';'810313';'812952';'813346';
          '822794'; '833557'; '841812'; '842426';'867237';'876501';'879873'];
      
%% 2. load class 1
num_cols = size(train_case, 1);
fileName = 'middle_linking_structure.m3d';
featureMatrix = [];
disp(['Loading features from all linking structures...']);
for i = 1 : size(train_case, 1)
    s_rep_path = fullfile(data_path, train_case(i, :));
    %exclude the link_to feature by setting the second para as 1
    featureVector = loadLinkingStructure(s_rep_path, fileName);    
    featureMatrix = [featureMatrix, featureVector];
end

%% 3. load class 2
fileName = 'diseased_linking_2.m3d';
for i = 1 : size(train_case, 1)
    s_rep_path = fullfile(data_path, train_case(i, :));
    featureVector = loadLinkingStructure(s_rep_path, fileName);    
    featureMatrix = [featureMatrix, featureVector];
end

%% 4. labels
y = ones(1, num_cols);
y(1:num_cols) = -1;
