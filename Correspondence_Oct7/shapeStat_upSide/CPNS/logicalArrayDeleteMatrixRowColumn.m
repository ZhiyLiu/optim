% x=rand(10);
% x([1,4,6,8],:) = 0;
% x([3,7],[2])=0;%% test and found: one zero in a row will be deleted????
% 
% disp(x);
% 
% % all(x,2) is a column vector, 1: if all elements in a row is zero, 0: if not.
% % index is a logical array, which can be used in logical indexing. Contain reverse side of all(x,2).  
% index = ~all(x,2); 
% disp(index);
% %% Here the logical array is different from index = [1;0;0;1;0;1;0;1;0;0];
% % Most arithmetic operations remove the logicalness from an array. For example, adding zero to a logical array removes its logical characteristic. 
% % A = +A is the easiest way to convert a logical array, A, to a numeric double array. Logical arrays are also created by the relational 
% % operators (==,<,>,~, etc.) and functions like any, all, isnan, isinf, and isfinite.
% 
% % convert the logical array into row vector, used to delete correspondence columns.
% index_c = ~(all(x,2)');
% disp(index_c);
% 
% % delete rows with all zero.
% x(index,:) = [];
% 
% % delete correspondence column.
% x(:,index_c)=[];
% 
% disp(x);



%% test
x=rand(10,7);
x([1,4,6,8],:) = 0;
x([3,7],[2])=0;%% test and found: one zero in a row will be deleted????

disp(x);

 index = all(x==0,2); 

% delete rows with all zero.
x(index,:) = [];

% x( all( ~any( x), 2 ), : ) = []; % removes all rows with all zero
% x( :, all( ~any( x ), 1 ) ) = []; % and columns

disp(index);



disp(x);