% This function plot the lines of area, volume, length change for up or down side while
% changing the interpoltaion level.
% num: is the subfolder number in each folder (e.g. if capture info from interpolation 0-8, the num is 9).
% dataDir: the folder where the different levels' .txt files stored.
% Command: plotInterpolationLevel('/NIRAL/work/ltu/WorkSpace/pictures/Fig-10-Interpolationlevel/up_or_down/up');

function [  ] = plotInterpolationLevel(dataDir, num)

%% Read in the datas (store in many .txt files) in the first subfolder, to get the size.
% Connect string to get the correct path
path = [dataDir, '/level_0'];
disp(path);

% Loop each file in this folder.
filename = strcat(path, '/boundaryHorizonal.txt');
disp(filename);
bh = load(filename);
[bh_r, bh_c] = size(bh);

filename = strcat(path, '/boundaryVertical.txt');
bv = load(filename);
[bv_r, bv_c] = size(bv);

filename = strcat(path, '/boundaryQuadsArea.txt');
bqa = load(filename);
[bqa_r, bqa_c] = size(bqa);

% Matrix used to store the 7 type of features.
boundaryHorizonal = zeros( num, bh_c) ;
boundaryVertical = zeros( num, bv_c) ;
boundaryQuadsArea = zeros( num, bqa_c) ;
quadsVolume = zeros( num, bqa_c) ;
skeletalHorizonal = zeros( num, bh_c) ;
skeletalVertical = zeros( num, bv_c) ;
skeletalQuadsArea = zeros( num, bqa_c) ;

%% Read in the datas (store in many .txt files) into one matrix
X = [0,1,2,3,5,8];
disp(X(1,1));
for i = 1:6,%1:num,
    j = X(1,i);
    % Connect string to get the correct path
    path = [dataDir, '/level_', int2str(j)];
    disp(path);

    % Loop each file in this folder.
    filename = strcat(path, '/boundaryHorizonal.txt');
    disp(filename);
    boundaryHorizonal(i, :) = load(filename);

    filename = strcat(path, '/boundaryVertical.txt');
    disp(filename);
    boundaryVertical(i, :) = load(filename);

    filename = strcat(path, '/boundaryQuadsArea.txt');
    disp(filename);
    boundaryQuadsArea(i, :) = load(filename);

    filename = strcat(path, '/quadsVolume.txt');
    disp(filename);
    quadsVolume(i, :) = load(filename);

    filename = strcat(path, '/skeletalHorizonal.txt');
    disp(filename);
    skeletalHorizonal(i, :) = load(filename);

    filename = strcat(path, '/skeletalVertical.txt');
    disp(filename);
    skeletalVertical(i, :) = load(filename);

    filename = strcat(path, '/skeletalQuadsArea.txt');
    disp(filename);
    skeletalQuadsArea(i, :) = load(filename);
end

% disp('boundaryHorizonal is: ');
% disp(boundaryHorizonal);
% disp('boundaryVertical is: ');
% disp(boundaryVertical);
% disp('boundaryQuadsArea is: ');
% disp(boundaryQuadsArea);
% disp('quadsVolume is: ');
% disp(quadsVolume);
% disp('skeletalHorizonal is: ');
% disp(skeletalHorizonal);
% disp('skeletalVertical is: ');
% disp(skeletalVertical);
% disp('skeletalQuadsArea is: ');
% disp(skeletalQuadsArea);

%% Save figure
figure(1);
plotAndSaveFig(boundaryQuadsArea, 'boundaryQuadsArea');
figure(2);
plotAndSaveFig(quadsVolume, 'quadsVolume');


end