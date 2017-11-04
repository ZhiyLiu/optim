% featurefilepath: /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepareavolume_down/
% path: /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepareavolume_down/*.txt 
% num: the eigenvalue number we want to keep to calculate the entropy.

function totalentropy = regularityEntropy(featurefilepath, num)

% Return the folder listing, restricted to files with a .txt extension, to the variable path:
path = strcat(featurefilepath, '*.txt');
files = dir(path);
% disp('Message from regularityEntropy.m; files in: ');
% disp(files);

% store the sum of each srep entropy
totalentropy = 0;

% Loop each file
for file = files'
%    disp('go into this file...');
%    disp(file);
   
   % Concatenating path and file name to get a full path of each feature file.
   filename = strcat(featurefilepath, file.name);
%    disp(filename);
   
   entropy = calculateEntropy(filename, num);
%    disp(entropy);
   
   totalentropy = totalentropy + entropy;
end
   
% disp(totalentropy);