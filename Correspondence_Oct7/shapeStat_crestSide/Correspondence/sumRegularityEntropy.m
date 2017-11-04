% Two return values: upentropy and downentropy

function [upentropy,downentropy] = sumRegularityEntropy()
%TOTALREGULARITYENTROPY Summary of this function goes here
%   Detailed explanation goes here


% For the down spokes:
% Calculate the srep area and volume entropy:
downareavolume = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepareavolume_down/',3);

% Calculate the srep horizonal edges entropy:
downhorizonal = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srephorizonal_down/',2);

% Calculate the srep vertical edges entropy:
downvertical = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepvertical_down/',2);

downentropy = downareavolume + downhorizonal + downvertical;
disp('The down entropy is: ');
disp(downentropy);

% For the top spokes:
% Calculate the srep area and volume entropy:
upareavolume = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepareavolume_top/',3);

% Calculate the srep horizonal edges entropy:
uphorizonal = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srephorizonal_top/',2);

% Calculate the srep vertical edges entropy:
upvertical = regularityEntropy('/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/output_feature/srepvertical_top/',2);

upentropy = upareavolume + uphorizonal + upvertical;
disp('The up entropy is: ');
disp(upentropy);

end

