%% ...... Test only ......
%     % Load a 3x10x5 matrix
%     load('/NIRAL/work/ltu/WorkSpace/shape/testData.mat', 'tA'); 
%      
%     [tm,tk,tn] = size(tA);
%     
%     % Create a weight matrix
%     W = zeros(1,tk,tn);    
%     
%     % Set some value to tW
%     for i=1:tn 
%         for j=1:tk
%             W(1,j,i) = randi(5);
%         end
%     end
%     
%     % Save A into a .mat file.
%     path = strcat('/NIRAL/work/ltu/WorkSpace/shape/', 'W.mat');
%     save(path,'W');
%     
%     load('/NIRAL/work/ltu/WorkSpace/shape/W.mat', 'tW'); 
%     
%     srepOutputData_GPA = GPA(tA, tW);
close all;
mcolor = [1 0 0];
markerstr = ['o'];
titlestr = {{'CPNS scores'}};
legendstr = {{''}};
[U S V] = svd(ZComp);
drawScatPlot(V, mcolor, legendstr, markerstr, titlestr);