function [x,a] = show_aitken(m,n,x0)
%
%  DESCRIPTION :
%    Applies first Newton's method to (x-1)^m,
%    using n steps starting at x0, followed by
%    Aitken acceleration method on the sequence. 
%
%  ON ENTRY :
%    m        multiplicity of 1 in (x-1)^m;
%    n        number of new approximations;
%    x0       starting value for the iteration.
%
%  ON RETURN :
%    x        sequence of approximation to 1,
%             as vector of range 1..n+1;
%    a        result of Aitken acceleration,
%             as vector of range 1..n-1.
%


% x = zeros(1,n+1);
% x(1) = x0;
% for k=1:n 
%   deltax = (x(k)-1)^m/(m*(x(k)-1)^(m-1));
%   x(k+1) = x(k) - deltax;
% end;

%x=[4.000000000000000e+000, 3.316624790355400e+000, 3.103747667048789e+000, 3.034385495301738e+000,3.011440019426500e+000,3.003810919291192e+000,3.001270037597814e+000,3.000423315999866e+000,3.000141102014992e+000,3.000047033636304e+000];
x=[0, 0.0737254, 0.0548481];
a = aitken(x);