% Calculate the geometry entropy of up and down spokes. 
% Input:
%   dataDir: folder of the m3d files you want to compute (e.g.
%   '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/aftermovewithdeltauv/'),
%   under this folder there
%   optionPNS: specify how CPNS should operate; 0 use SMALL circles or BIG
%   circles depending upon analysis; 1 always use SMALL circles; 2 always
%   use BIG circles.
%
% Output:
%   upGE: up geometry entropy
%   downGE: down geometry entropy
 
function [ upGE, downGE ] = calculateGeometryEntropy( dataDir, optionPNS )

upDataDir = strcat(dataDir, 'up/');
disp(strcat('The files for computing up entropy is in: ', upDataDir));
% for up spokes
upGE = calculateCPNS_eigenvalue(upDataDir, optionPNS, 0);
disp('The up geometry entropy is: ');
disp(upGE);


downDataDir = strcat(dataDir, 'down/');
disp(strcat('The files for computing down entropy is in: ', downDataDir));
% for down spokes
downGE = calculateCPNS_eigenvalue(downDataDir, optionPNS, 1);
disp('The down geometry entropy is: ');
disp(downGE);

end

