function Logpx = LogNPd(x)
% LOGNP Riemannian log map at North pole of S^k
%       LogNP(x) returns k x n matrix where each column is a point on tangent
%       space at north pole and the input x is (k+1) x n matrix where each column 
%       is a point on a sphere.
%
%
%   See also ExpNPd.

% Last updated Oct 20, 2009
% Sungkyu Jung



% Liyun Tu Dec 2, 2013
% Logpx will be a complex number if x(end,:).^2 is bigger than one, because
% minus number use sqrt will be a complex number.
% acos(-1.0009) is complex number.

% x is:
%    0.0266    0.0261    0.0269    0.0239    0.0279    0.0257    0.0296    0.0245    0.0254    0.0250
%    0.0178   -0.0173    0.0175    0.0470    0.0020   -0.0203   -0.0223    0.0363    0.0253    0.0065
%   -0.9989   -0.9988   -1.0009   -0.9993   -1.0001   -1.0010   -0.9998   -1.0006   -1.0000   -0.9996

% x(end,:) is:
%   -0.9989   -0.9988   -1.0009   -0.9993   -1.0001   -1.0010   -0.9998   -1.0006   -1.0000   -0.9996

% acos(x(end,:)) is:
%  3.0957     3.0935    3.1416 - 0.0423i   3.1048    3.1416 - 0.0117i   3.1416 - 0.0439i   3.1206   3.1416 - 0.0360i   3.1396   3.1142  

% x(end,:).^2 is:
%    0.9979    0.9977    1.0018    0.9987    1.0001    1.0019    0.9996    1.0013    1.0000    0.9993

% 1-x(end,:).^2
%    0.0021    0.0023   -0.0018    0.0013   -0.0001   -0.0019    0.0004   -0.0013    0.0000    0.0007

% Message from LogNPd.m;
% x(1:(end-1),:) is:
%    0.0266    0.0261    0.0269    0.0239    0.0279    0.0257    0.0296    0.0245    0.0254    0.0250
%    0.0178   -0.0173    0.0175    0.0470    0.0020   -0.0203   -0.0223    0.0363    0.0253    0.0065







% disp('x is: ');
% disp(x);
[d n] = size(x);
%scale = acos(x(end,:))./sqrt(abs(1-x(end,:).^2));
% disp(acos(x(end,:)));
scale = acos(x(end,:))./sqrt(1-x(end,:).^2);
% 
% disp('1-x(end,:).^2');
% disp(1-x(end,:).^2);
% disp('abs(1-x(end,:).^2)');
% disp(abs(1-x(end,:).^2));

scale(isnan(scale)) =1;
% disp('Message from LogNPd.m; ');
% disp(scale);
%disp(repmat(scale,d-1,1));
Logpx = repmat(scale,d-1,1).*x(1:(end-1),:);


