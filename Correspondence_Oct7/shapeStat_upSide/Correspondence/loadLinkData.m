%% 1. collect all links from 10 cases, each case relates to 1 column
% s_reps = ['/playpen/workspace/demo_multiobject/180215';
%     '/playpen/workspace/demo_multiobject/223443';
%     '/playpen/workspace/demo_multiobject/225196';
%     '/playpen/workspace/demo_multiobject/402809';
%     '/playpen/workspace/demo_multiobject/549967';
%     '/playpen/workspace/demo_multiobject/622437';
%     '/playpen/workspace/demo_multiobject/721965';
%     '/playpen/workspace/demo_multiobject/730231';
%     '/playpen/workspace/demo_multiobject/802785';
%     '/playpen/workspace/demo_multiobject/867237'];
% 
% legend_str = {{'180215','223443','225196','402809','549967','622437','721965','730231','802785','867237',...
%     '180215_diseased','223443_diseased','225196_diseased','402809_diseased',...
%     '549967_diseased','622437_diseased','721965_diseased','730231_diseased',...
%     '802785_diseased','867237_diseased'}};

data_path = '/playpen/workspace/proj_multiobject/';
s_reps = ['107524'; '141335'; '145680'; '180215';'223443';'225196';'236316';'252848';'315149';'323623';
          '339961'; '359973'; '360222'; '365557';'402809';'429591';'482642';'490181';'513330';'549587';
          '549967'; '589344'; '618672'; '622437';'631202';'638999';'641078';'660372';'691501';'699209';
          '713824'; '721965'; '730231'; '755016';'792210';'793001';'802785';'810313';'812952';'813346';
          '889945'; '895096'; '910593'; '915717';'931663';'938659';'950194';'961753';'967950';'974849';
          '822794'; '833557'; '841812'; '842426';'867237';'876501';'879873'];
% 12-dimensional(purlz+linkTo) links and 68 links per object
% each column consists of feature vector per configuration(3 objects)

num_rows = 12 * 68 * 3; % == 2448
num_cols = size(s_reps, 1);
fileName = 'middle_linking_structure.m3d';
shapeFeature = [];
linkFeature = [];
disp(['Loading features from all linking structures...']);
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(data_path, s_reps(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    shapeFeature = [shapeFeature, shape];
    linkFeature = [linkFeature, link];
end

% load 1st disease group
% fileName = 'diseased_linking_structure.m3d';
% for i = 1 : size(s_reps, 1)
%     s_rep_path = fullfile(s_reps(i, :),fileName);
%     [shape link] = loadJIVE(s_rep_path);
%     shapeFeature = [shapeFeature, shape];
%     linkFeature = [linkFeature, link];
% end


% load 2nd 'diseased' linkinking structure
fileName = 'diseased_linking_2.m3d';
for i = 1 : size(s_reps, 1)
    s_rep_path = fullfile(data_path, s_reps(i, :),fileName);
    [shape link] = loadJIVE(s_rep_path);
    shapeFeature = [shapeFeature, shape];
    linkFeature = [linkFeature, link];
end

%% Exclude sensitive features result from observations on loading plots
% idSensitiveFeature = [37]; % reversely ordered
% sizePerSpoke = 7;
% sizePerLink = 5;
% numSamples = size(shapeFeature, 1);
% for i = 1 : length(idSensitiveFeature)
%     id = idSensitiveFeature(i);
%     shapeFeature((id -1) * sizePerSpoke+1: id*sizePerSpoke, :) = [];
%     linkFeature((id -1) * sizePerLink+1: id*sizePerLink, :) = [];
% end