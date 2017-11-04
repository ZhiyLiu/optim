function [sym] = TubeSymPrimitive( varargin )
%
% [prim] = TubeSymPrimitive()
%
% Constructs a new default tube primitive object
%
% [prim] = TubeSymPrimitive(prim2)
%
% Copy constructor
%
% [prim] = TubeSymPrimitive( pos, r, elongation, U0 ,halfConeAngle, baseAtom, dr [, inTangentSpace])
%
% Constructs a symmetric-space tube primitive by using the values supplied as parameters.
%

switch(nargin)
	case 0		% default constructor
		s	= struct( ...
			'pos', [ 0.0; 0.0; 0.0 ], ...	% the position of the hub in 3 space
			'r', 1.0, ...					% the length of the spokes
			'elongation', 1.0, ...			% the elongation if an end primitive.
			'U0', [ 1.0; 0.0; 0.0 ], ...	% the direction of the tangent vector
			'hca', pi/2.0, ...				% the half cone angle
			'baseAtom', false, ...			% Whether base atom or not
			'dr', repmat( 0.0, [8,1] ), ...	% Deviations of the 8 spokes.
			'inTangentSpace', false );		% whether we are in tangent space or not
		sym	= class( s, 'TubeSymPrimitive' );
	case 1		% copy constructor
		if( isa(varargin{1}, 'TubeSymPrimitive') )
			sym	= varargin{1};
		else
			error([ 'Argument 1 should be of type TubeSymPrimitive, it is instead of type ' class(varargin{1}) ]);
		end
	case 7
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'U0', varargin{4}, ...			% the direction of the tangent vector
			'hca', varargin{5}, ...			% the half cone angle
			'baseAtom', varargin{6}, ...	% whether base atom or not
			'dr', varargin{7}, ...	% Deviations of the 8 spokes.
			'inTangentSpace', false );		% whether we are in tangent space or not
		sym	= class(s, 'TubeSymPrimitive');
	case 8
		s	= struct( ...
			'pos', varargin{1}, ...			% the position of the atom in 3-space
			'r', varargin{2}, ...			% the spoke length
			'elongation', varargin{3}, ...	% the elongation if an end atom
			'U0', varargin{4}, ...			% the direction of the tangent vector
			'hca', varargin{5}, ...			% the half cone angle
			'baseAtom', varargin{6}, ...	% whether base atom or not
			'dr', varargin{7}, ...			% Deviations of the 8 spokes.
			'inTangentSpace', varargin{8});	% whether we are in tangent space or not
		sym	= class(s, 'TubeSymPrimitive');
	otherwise
		error('Incorrect constructor call, see help');
end

return;
