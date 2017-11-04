%% u direction as axis
% Rotation angle
theta = -pi/12;

%ux, uy, uz is the direction of the rotation axis.
ux = 1;
uy = 0;
uz = 0;

% Calculate the rotation matrix
costheta = cos(theta);
v1 = costheta + ux*ux*(1-costheta);
v2 = -(ux*sin(theta))
R = [v1,0,0;0,costheta,v2;0,-v2,costheta];

% The spoke direction need to be rotated. 3*3 srep, primitive[1,1]'s up
% spoke.
U = [0,0,-1];
%U = [ 0,    -2.588190451025207e-01,    -9.659258262890683e-01];

result = R*U'


%% v direction as axis
% % Rotation angle
% theta = pi/7;
% 
% %ux, uy, uz is the direction of the rotation axis.
% ux = 0;
% uy = 1;
% uz = 0;
% 
% % Calculate the rotation matrix
% costheta = cos(theta);
% v1 = costheta + uy*uy*(1-costheta);
% v2 = -(uy*sin(theta))
% R = [costheta,0,sin(theta);0,v1,0;v2,0,costheta];
% 
% % The spoke direction need to be rotated. 3*3 srep, primitive[1,1]'s up
% % spoke.
% U = [0,0,-1];
% %U = [2.2204460492503131e-16, 6.1232339957367623e-17, -1];
% 
% result = R*U'