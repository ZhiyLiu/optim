% Constructor of class 'QuadPrimitive'

function quadPrim = QuadPrimitive(varargin)

% quadPrim = QuadPrimitive(pos: 3x1 vector, r, elongation, q: 4x1(or 3x1) vector, theta (in radian),
%                          inTangentSpace: boolean)
% NOTE : q = [qw; qx; qy; qz] qw = cos(theta/2), [qx, qy, qz] = sin(theta/2)*(axis of rotation)

% 	struct( ...
% 			'pos', [ 0.0; 0.0; 0.0 ], ...	% the position of the atom in 3-space
% 			'r', 1.0, ...					% the spoke length
% 			'elongation', 1.0, ...			% the elongation if an end atom
% 			'q', [ 1.0; 0.0; 0.0; 0.0], ...	% a quaternion representing the orientation
% 			'theta', pi/2.0, ...			% the half angle between the spokes
% 			'inTangentSpace', false );		% whether we are in tangent space or not (always false for tubes)

switch nargin
    case 0  % default (quadPrimitive in manifold)
        quadPrim = struct('pos', [0.0; 0.0; 0.0], 'r', 1.0, 'elongation', 1.0,...
            'q', [1.0; 0.0; 0.0; 0.0], 'theta', pi/2.0, 'inTangentSpace', false);
        quadPrim = class(quadPrim, 'QuadPrimitive');
    case 1
        if (isa(varargin{1}, 'QuadPrimitive'))
            quadPrim = varargin{1};
        else
            error('Wrong argument type');
        end
    case 2
        quadPrim =	struct( ...
			'pos', [ 0.0; 0.0; 0.0 ], ...	% the position of the atom in 3-space
			'r', 1.0, ...					% the spoke length
			'elongation', 1.0, ...			% the elongation if an end atom
			'q', [ 1.0; 0.0; 0.0; 0.0], ...	% a quaternion representing the orientation
			'theta', pi/2.0, ...			% the half angle between the spokes
			'inTangentSpace', false );		% whether we are in tangent space or not (always false for tubes)
        
        primKeyStr	= varargin{1};
        modelArray	= varargin{2};
        quadPrim.pos = [ ...
            findVal([primKeyStr 'x'], modelArray, ' %f'); ...
            findVal([primKeyStr 'y'], modelArray, ' %f'); ...
            findVal([primKeyStr 'z'], modelArray, ' %f') ];
        quadPrim.r  = findVal([primKeyStr 'r'], modelArray, ' %f');
        type		= findVal([primKeyStr 'type'], modelArray, ' %s' );
        if( strcmp( type, 'EndPrimitive' ))
            quadPrim.elongation	= findVal([primKeyStr 'elongation'], modelArray, ' %f' );
        elseif( ~strcmp( type, 'StandardPrimitive' ))
            %error( ['Unknown primitive type' type ]);
            quadPrim.elongation	= 1.0;
        end
        q		= [ ...
            findVal([primKeyStr 'qw'], modelArray, ' %f'); ...
            findVal([primKeyStr 'qx'], modelArray, ' %f'); ...
            findVal([primKeyStr 'qy'], modelArray, ' %f'); ...
            findVal([primKeyStr 'qz'], modelArray, ' %f') ];
        %% Make sure qw>0
        if (q(1) < 0)
            q(1) = -q(1);
            q(2) = -q(2);
            q(3) = -q(3);
            q(4) = -q(4);
        end
        quadPrim.q  = q;
        theta	= findVal([primKeyStr 'theta'], modelArray, ' %f');
        % Angle is saved in degrees, convert to radians.
        quadPrim.theta	= theta * pi/180.0;
        quadPrim = class(quadPrim, 'QuadPrimitive');
    case 5
        if (length(varargin{4}) == 3)
            quadPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
                'q', varargin{4}, 'theta', varargin{5}, 'inTangentSpace', true);
        elseif (length(varargin{4}) == 4)
            q =  varargin{4};
            %% Make sure qw>0
            if (q(1) < 0)
                q(1) = -q(1);
                q(2) = -q(2);
                q(3) = -q(3);
                q(4) = -q(4);
            end
            quadPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
                'q', q, 'theta', varargin{5}, 'inTangentSpace', false);
        else
            error('Wrong argument');
        end
        quadPrim = class(quadPrim, 'QuadPrimitive');
    case 6
        inTangentSpace = varargin{6};
        if (~inTangentSpace)
            q =  varargin{4};
            %% Make sure qw>0
            if (q(1) < 0)
                q(1) = -q(1);
                q(2) = -q(2);
                q(3) = -q(3);
                q(4) = -q(4);
            end
        end     
        quadPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
            'q', q, 'theta', varargin{5}, 'inTangentSpace', varargin{6});
        quadPrim = class(quadPrim, 'QuadPrimitive');
    otherwise
        error('Wrong number of input arguments')
end



