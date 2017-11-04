% =====================================================================
% function [COG, T, S, R, Q] = readSim(dirPath)
%
%   Read 'estimated_similarity.data' file that prediction program 
%   spits out. The file contains the similarity transform in the 
%   following format(without letters and ())
%
%   center_of_rotation (X, Y, Z) : COG
%   translation (X, Y, Z) : T
%   scale (S) : S
%   rotation_matrix (3 x 3) : R
%   quarternion (W, X, Y, Z) : Q
% ---------------------------------------------------------------------
function [COG, T, S, R, Q] = readSim(dirPath)

fileName = 'estimated_similarity.dat';
filePath = fullfile(dirPath, fileName);
%format long e;
%I think, although the display of the matlab looks like it rounds up the
%numbers that are read in, the actual numbers stored in the variable
%retain the precision.
data = dlmread(filePath);
COG = data(1,1:3)';
T = data(2,1:3)';
S = data(3,1);
R = data(4:6, 1:3)';
Q = data(7, :)';

delete(fileName);





