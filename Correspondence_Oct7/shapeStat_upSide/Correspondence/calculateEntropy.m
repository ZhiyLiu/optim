% num is the eigenvalue number used.
% filename is the feature file store the feature matrix

function entropy = calculateEntropy(filename, num)
%filename = '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/popfeaturematrix_top.txt';
[A,delimiterOut]=importdata(filename);
 disp('The input A matrix is: ');
 disp(A);

% A is a matrix with row as feature(dimesion), column as a case(sample).
[m,n] = size(A);

% mean for each dimension(one row is one dimension), its a m by 1 matrix.
C = mean(A,2);

% subtract off the row mean from A
M = repmat(C,1,n); %repmat(C,1,n) means tile 1*n patch of C.
centereddata = A - M;
disp('The centereddata matrix is: ');
disp(centereddata);

%% add fixed noise
theta = 0;%0.001;
rng('default');
du = rand(m, n);
% disp('The du is: ');
% disp(du);
centereddata = centereddata + theta*du;
%   disp('The centereddata matrix after advising is: ');
%   disp(centereddata);

% calculate the covariance matrix
trancentereddata = transpose(centereddata);
covariance_1 = 1/(n-1) * centereddata * trancentereddata;
% make the covariance matrix symmetric.
covariance = (covariance_1 + transpose(covariance_1))/2;
 disp('The covariance matrix is: ');
 disp(covariance);
%  dlmwrite('/NIRAL/work/ltu/WorkSpace/covMatrix.txt',covariance,' ');

%matlab presice cause the covariance is not a symmetric matrix.
%%temp = covariance - covariance_1;%the result is not 0, with very small numbers due to the Data Accuracy;
%deside whether the covariance matrix is a symmetric matrix.
%%p = transpose(covariance);
%%pp=covariance-p;
%%o = any(pp(:)); %0: symmetric; 1: not symmetric.

% test 
%X = cov(trancentereddata);%matlab cov function gives the same result as covariance.

% Convert covariance matrix to correlation matrix.
correlation_1 = corrcov(covariance);
ttt = correlation_1 - transpose(correlation_1);
correlation = (correlation_1 + transpose(correlation_1))/2;
 disp('The correlation matrix is: ');
 disp(correlation);


% find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V)
[PC, V] = eig(correlation);
%  disp('V is: ');
%  disp(V);
% extract diagonal of matrix (eigenvalues) as vector
V = diag(V);

% sort the variances in decreasing order
% rindices is the original index of element in V, use this we can get the
% element from V in decreasing order.
[junk, rindices] = sort(-1*V);
eigenvalue = V(rindices); 
%%eigenvector = PC(:,rindices);%for all row, order element using column index rindices.

%gg = [0.677873398528012,0.735178655544408;-0.735178655544408,0.677873398528012];
%%gg = transpose(eigenvector);
% project the original data set
%%finalData = gg * centereddata;

% get back the original data set 
%%RowOriginalData = eigenvector * finalData + M;

%jj = eigenvalue(1,1);
%jjj = eigenvalue(2,1);

%%mmm = log(-2);
%  disp('eigenvalue Reg');
%  disp(eigenvalue);
%  dlmwrite('/NIRAL/work/ltu/WorkSpace/eigenvalue.txt',eigenvalue,' ');
% 
% for i = 1:num    %only use num eigenvalue. 
%     entropy = sum(log(eigenvalue(i,1)));
% end

    threshold = eigenvalue/sum(eigenvalue); 

    eigenvalue(threshold<0.01) = [];
    
    entropy = sum(log(eigenvalue));










