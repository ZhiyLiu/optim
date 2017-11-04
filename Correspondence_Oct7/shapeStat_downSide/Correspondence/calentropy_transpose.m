% Reference: http://www.cs.stevens.edu/~mordohai/classes/cs559_s09/PCA_in_MATLAB.pdf
% This first version follows Section 5 (of above tutorial) by examining the covariance of the data set.


% PCA1: Perform PCA using covariance.

% load the training set features from a txt file into a matrix A. 
% Column is case numbers, row is features.
filename = '/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/test_PCA1.txt';
[A,delimiterOut]=importdata(filename);

% B is a matrix with row as a case(sample), each column is a feature(dimesion).
B = transpose(A);

[m,n] = size(B);

% mean for each dimension(column), its a 1 by n matrix.
C = mean(B,1);

% subtract off the column mean from B
M = repmat(C,m,1); %repmat(C,m,1) means tile m*1 patch of C.
centereddata = B - M;  

% calculate the covariance matrix
trancentereddata = transpose(centereddata);
covariance = 1/(m-1) * trancentereddata * centereddata;

% test 
%X = cov(B);%matlab cov function gives the same result as covariance.

% find the eigenvectors(PC, column vector, each column is a PC) and eigenvalues(diagonal of V)
[PC, V] = eig(covariance);

% extract diagonal of matrix (eigenvalues) as vector
V = diag(V);

% sort the variances in decreasing order
% rindices is the original index of element in V, use this we can get the
% element from V in decreasing order.
[junk, rindices] = sort(-1*V);
eigenvalue = V(rindices); 
eigenvector = PC(:,rindices);%for all row, order element using column index rindices.


%gg = [0.677873398528012,0.735178655544408;-0.735178655544408,0.677873398528012];
gg = transpose(eigenvector);
% project the original data set
finalData = gg * trancentereddata; 

RowOriginalData = transpose(gg) * finalData + transpose(M);