function [ output_args ] = normalDistributionMergedudv( input_args )

% variables for du
A = load('/NIRAL/work/ltu/WorkSpace/Test/Apr22/norDisNums-100-du.txt');
[m,n] = size(A);

% variables for dv
B = load('/NIRAL/work/ltu/WorkSpace/Test/Apr22/norDisNums-100-dv.txt');

% Merge A and B to get du, dv, du, dv, ...
k = 2*n
C = zeros(m, k);
for i = 1:n,
    C(m, 2*(i-1)+1) = A(m,i); % set value to C matrix entries: 1, 3, 5, 7, ..., 2*(i-1)+1
    C(m, 2*i) = B(m, i);
end

% write to file
dlmwrite('/NIRAL/work/ltu/WorkSpace/Test/Apr22/merge.txt', C, ' ');


end

