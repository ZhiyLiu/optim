function [q] = quatFromFrame(tangent, normal)
%
% [q] = quatFromFrame( tangent, normal )
%
% This functions constructs a unit quaternion representing
% a rotation of the standard frame into the frame represented
% by [tangent, normal, tangent x normal ]
% tangent, normal should be unit vectors, no check is enforced.
%

%
% NOTE: IF you find any fixes in this code, then change the code in Pablo
% in Quat::buildFromFrame also.
%

R(1:3)	= tangent;
R(4:6)	= normal;
R(7:9)	= cross(tangent,normal);

T = R(1)+R(5)+R(9);
% @Jeong
% NOTE: T+1 == 0 when rotation by 180 degree == flip 
% In converting rotation matrix to quaternion,
% we need to handle this case specially.   
if ( (1+T) > 1e-14 )
	cos_phi2 = sqrt(T+1)/2; % (==) cos_phi2 = cos( acos( ((R(1)+R(5)+R(9))-1)/2 ) / 2 );
	rot_axis = [R(6)-R(8); R(7)-R(3); R(2)-R(4)] ./ repmat(4*cos_phi2,3,1);            
else  % Flip
	cos_phi2 = 0.0;
	if ( (1+R(1)) > 0)
		qx = sqrt(0.5*(1+R(1)));
		S =0.5/qx;
		qy = R(4)*S;
		qz = R(7)*S;   
	elseif ( (1+R(5)) > 0)
		qy = sqrt(0.5*(1+R(5)));
		S = 0.5/qy;
		qx = R(4)*S;
		qz = R(8)*S;
	elseif ( (1+R(9)) > 0)
		qz = sqrt(0.5*(1+R(9)));
		S = 0.5/qz;
		qx = R(7)*S;
		qy = R(8)*S;
	else
		error(' In quatFromFrame, this error should not happen.')
	end
	rot_axis = [qx; qy; qz];
end
q = [cos_phi2; rot_axis];
q_norm = sqrt(sum(q.^2));
q = q./q_norm;

