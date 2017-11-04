%% test1
% A = load('/NIRAL/work/ltu/WorkSpace/Test/Feb20_forpaper/test1/merge.txt');
% 
% B = zeros(1,100);
% 
% k=1;    
%     for j = 1 : 50,
%         for i = 1 : 2,            
%             B(1, k) = A(i,j);
%             k=k+1;
%         end
%     end    
% 
% dlmwrite('/NIRAL/work/ltu/WorkSpace/Test/Feb20_forpaper/test1/dudv.txt',B,' ');




%% test2
% A = load('/NIRAL/work/ltu/WorkSpace/Test/Feb20_forpaper/o.txt');
% 
% B = zeros(1,120);
%  
% 
% k=1;    
%     for j = 1 : 20,
%         for i = 1 : 6,            
%             B(1, k) = A(i,j);
%             k=k+1;
%         end
%     end    
% 
% dlmwrite('/NIRAL/work/ltu/WorkSpace/Test/Feb20_forpaper/dudv.txt',B,' ');


%% test-- May23-200Up
A = load('/NIRAL/work/ltu/WorkSpace/Test/May23-200Up/input/merge.txt');

B = zeros(1,400);

k=1;    
    for j = 1 : 200,
        for i = 1 : 2,            
            B(1, k) = A(i,j);
            k=k+1;
        end
    end    

dlmwrite('/NIRAL/work/ltu/WorkSpace/Test/May23-200Up/input/dudv.txt',B,' ');