diary output_aitken;
format long e
[x,a] = show_aitken(3,5,1.001);
% fprintf('original sequence ');
% x = x';
% fprintf('accelerated sequence ');
% a = a';
% fprintf('original errors ');
% e = x - ones(size(x));
% fprintf('errors after acceleration ');
% error = a - ones(size(a));
diary off;