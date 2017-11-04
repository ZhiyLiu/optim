function [ output_args ] = normalDistributionVars( input_args )
%NORMALDISTRIBUTIONVARS Summary of this function goes here
%   Detailed explanation goes here

% numSamples = 30;
% mu = 0;
% sigma = 1;
% samples = mu + sigma.*randn(numSamples, 1);
% %samples = randn(numSamples, 1);
% disp(samples);
% figure;hist(samples(:));

%% rand() is Uniformly distributed pseudorandom numbers. from 0 to 1.
% clear du
% tic
% du = rand(3, 10);
% disp(du);
% disp('The mean of the generated random vector is: ');
% disp(mean(mean(du)));
% toc

%% randn() is Normally distributed pseudorandom numbers. May over out 0 to 1.
% clear du
% tic
% du=randn(3, 10);
% disp(du);
% disp('The mean of the generated random vector is: ');
% disp(mean(mean(du)));
% toc


%% Generate normally distributed numbers between (0,1).
a=0+1*randn(1,3000); %generate a row vector with 1000 normally distributed number with mean=0, standard_deviation=1.
a=a((a>-0.5 & a<0.5));%select those between (-0.4,0.4)
a=a(1:28);%select the first 50 ones.
disp('The numbers generated is: ');
disp(a);
figure;hist(a(:));
dlmwrite('/NIRAL/work/ltu/WorkSpace/norDisNums-28_crest.txt',a,'delimiter',' ','precision', 2);% set precision, keep 16....


end

