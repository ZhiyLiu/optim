function [prim] = TubePrimitive(varargin);
%
% [prim] = TubePrimitive()
%
% Constructs a new default tube primitive object
%
% [prim] = TubePrimitive(prim2)
%
% Copy constructor
%
% [prim] = TubePrimitive( primKeyStr, modelArray )
%
% Constructs a tube primitive by reading from the cell array of strings,
% modelArray. primKeyStr is the key to look for.
%
% [prim] = TubePrimitive( pos, r, elongation, q, theta [, inTangentSpace [,baseAtom [, dr] ] ])
%
% Constructs a tube primitive by using the values supplied as parameters.
%

switch(nargin)
	case 0		% default constructor
		s	= struct( ...
			'pos', [ 0.0; 0.0; 0.0 ], ...	% the position of the atom in 3-space
			'r', 1.0, ...					% the spoke length
			'elongation', 1.0, ...			% the elongation if an end atom
			'q', [ 1.0; 0.0; 0.0; 0.0], ...	% a quaternion representing the orientation
			'theta', pi/2.0, ...			% the half angle between the spokes
			'baseAtom', false, ...			% Whether this atom is a base atom or not.
			'dr', repmat( 0.0, [8,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', false );		% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
	case 1		% copy constructor
		if( isa(varargin{1}, 'TubePrimitive') )
			prim	= varargin{1};
		else
			error([ 'Argument 1 should be of type TubePrimitive, it is instead of type ' class(varargin{1}) ]);
		end
	case 3
		primKeyStr	= varargin{1};
		modelArray	= varargin{2};
        numberOfSpokes = varargin{3};
		s	= struct( ...
			'pos', [ 0.0; 0.0; 0.0 ], ...	% the position of the atom in 3-space
			'r', 1.0, ...					% the spoke length
			'elongation', 1.0, ...			% the elongation if an end atom
			'q', [ 1.0; 0.0; 0.0; 0.0], ...	% a quaternion representing the orientation
			'theta', pi/2.0, ...			% the half angle between the spokes
			'baseAtom', false, ...			% Whether this atom is a base atom or not.
			'dr', repmat( 0.0, [numberOfSpokes,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', false );		% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
		prim.pos	= [ ...
			findVal([primKeyStr 'x'], modelArray, ' %f'); ...
			findVal([primKeyStr 'y'], modelArray, ' %f'); ...
			findVal([primKeyStr 'z'], modelArray, ' %f') ];
		prim.r		= findVal([primKeyStr 'r'], modelArray, ' %f');
		type		= findVal([primKeyStr 'type'], modelArray, ' %s' );
		if( strcmp( type, 'EndPrimitive' ))
			prim.elongation	= findVal([primKeyStr 'elongation'], modelArray, ' %f' );
		elseif( ~strcmp( type, 'StandardPrimitive' ))
			error( ['Unknown primitive type' type ]);
		end
		prim.q		= [ ...
			findVal([primKeyStr 'qw'], modelArray, ' %f'); ...
			findVal([primKeyStr 'qx'], modelArray, ' %f'); ...
			findVal([primKeyStr 'qy'], modelArray, ' %f'); ...
			findVal([primKeyStr 'qz'], modelArray, ' %f') ];
		prim.theta	= findVal([primKeyStr 'theta'], modelArray, ' %f');
		% Angle is saved in degrees, convert to radians.
		prim.theta	= prim.theta * pi/180.0;
		prim.baseAtom	= findVal([primKeyStr 'baseAtom'], modelArray, ' %f');
		for i=1:numberOfSpokes
			prim.dr(i)	= findVal([primKeyStr 'dr[' num2str(i-1) ']' ], modelArray, ' %f');
		end
	case 5
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'q', varargin{4}, ...			% a quaternion representing the orientation
			'theta', varargin{5}, ...		% the half angle between the spokes
			'baseAtom', false, ...			% Whether this atom is a base atom or not.
			'dr', repmat( 0.0, [8,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', false );		% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
	case 6
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'q', varargin{4}, ...			% a quaternion representing the orientation
			'theta', varargin{5}, ...		% the half angle between the spokes
			'baseAtom', false, ...			% Whether this atom is a base atom or not.
			'dr', repmat( 0.0, [8,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', varargin{6});	% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
	case 7
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'q', varargin{4}, ...			% a quaternion representing the orientation
			'theta', varargin{5}, ...		% the half angle between the spokes
			'baseAtom', varargin{7}, ...	% Whether this atom is a base atom or not.
			'dr', repmat( 0.0, [8,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', varargin{6});	% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
	case 8
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'q', varargin{4}, ...			% a quaternion representing the orientation
			'theta', varargin{5}, ...		% the half angle between the spokes
			'baseAtom', varargin{7}, ...	% Whether this atom is a base atom or not.
			'dr', varargin{8}, ...			% Deviations of the 8 spokes.
			'inTangentSpace', varargin{6});	% whether we are in tangent space or not (always false for tubes)
		prim	= class(s, 'TubePrimitive');
	otherwise
		error('Incorrect constructor call, see help');
end

return;
