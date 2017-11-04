function [R] = quatToRotationMatrix(q)
%
% [R] = quatToRotationMatrix(q)
%
% Given a quaternion q, this function returns the corresponding rotation matrix R
%

R	= zeros(3,3);

%
% For understanding this code, refer to
% http://mathworld.wolfram.com/EulerParameters.html
%

qij	= q*q';

for i = 1:3
	for j = 1:3
		d	= i==j;
        
        R(i,j) = d * (q(1).^2 - sum(q(2:end).^2)) + 2 * qij(i+1, j+1);
        for k = 1:3
            % the permutation symbol.
            if( i == j || j == k || k == i )
                p	= 0;
            elseif( mod(j-i,3) == 1 && mod(k-j,3) == 1 && mod(i-k,3) == 1 )
                p	= 1;
            else
                p	= -1;
            end

            R(i,j)	= R(i,j) + 2 * p * qij(1,k+1);
        end
	end
end

