% This objective funciton considered both the geometry entropy and
% regularity entropy.

function [ output_args ] = entropy( input_args )
%ENTROPY Summary of this function goes here
%   Detailed explanation goes here


% The sum of N srep's regularity entropy. 
[upentropy,downentropy] = sumRegularityEntropy();
disp(upentropy);
disp(downentropy);


% The geometric entropy. Also seperated into up and down parts.



end

