function [tube] = alignTube(tube, alignWith)
%
% [tube] = alignTube(tube, alignWith)
%
% tube is an array of TubePrimitives that you want to be phi-re-aligned.
% The alignment is obtained as the average of the normals of all the atoms
% contained in the array alignWith that are marked as baseAtoms.
%
% The index of the first atom in the alignWith array is used to flag the base
% atom in the tube array.
%

if( mod(prod(size(alignWith)), length(tube)) ~= 0 )
	error('The tubes to be aligned with must be of the same length as the first one');
end

if( length(tube) == 1 )
	% Quietly ignore 1 atom long tubes
	return;
end

% Find average normal.
ibaseAtoms	= find([alignWith.baseAtom] == 1);

%
% We could put more assertions in here, such as the base atom stays the same across the
% entire population .....
%

if( length(ibaseAtoms) == 0 )
	ibaseAtoms	= 1:length(tube):prod(size(alignWith));
	ibaseAtoms	= ibaseAtoms + ceil(length(tube)/2);
end

baseAtomsq	= [alignWith.q];
baseAtomsq	= baseAtomsq(:,ibaseAtoms);
n	= [];
for i = 1:size(ibaseAtoms)
	% rotate vector [0,1,0] (the normal to the medial axis) into
	% position.
	n	= [n QuatRotVec( baseAtomsq(:,i), [0; 1; 0] )];
end

meanN	= SphereMean(n, eye(length(n)/3) );

%
% Set base atom
%
prim	= tube((ibaseAtoms(1)-1)/size(alignWith,1)+1);
T		= QuatRotVec( prim.q, [1; 0; 0] );
% make meanN orthogonal to T.
meanN	= meanN - dot(meanN, T) * T;
prim.q	= quatFromFrame( T, meanN );
prim.baseAtom	= true;
tube(ibaseAtoms(1))	= prim;

%
% Orient all the other normals.
%
for side=0:1
	if side == 0
		range	= ibaseAtoms(1)-1:-1:1;
	else
		range	= ibaseAtoms(1)+1:1:length(tube);
	end
	prevPrim	= tube(ibaseAtoms(1));
	for i = range
		prim	= tube(i);
		% Get tangent vectors
		T		= QuatRotVec( prim.q, [1; 0; 0] );
		prevT	= QuatRotVec( prevPrim.q, [1; 0; 0] );
		
		dT	= T - prevT;
		normdT	= sqrt(sum(dT.*dT));
		if( normdT < 1.0e-14 )
			% Both tangents are in almost perfect alignment
			% Just copy over quaternion.
			prim.q	= prevPrim.q;
		else
			% Use dT direction as the vector along which to minimize distances
			% Find a normal to T to use
			prevNormal	= dT - dot(dT,prevT) * prevT;
			prevNormal	= prevNormal / sqrt(sum(prevNormal.*prevNormal));
			normal		= dT - dot(dT,T)*T;
			normal		= normal / sqrt(sum(normal.*normal));
			%
			% Derivation:
			% pT = S [1 0 0]'		T = Q [1 0 0]'
			% pN = S [1 0 0]'		N = Q [0 1 0]'
			% S is prevPrim.q;
			% p means previous (one with fixed and known S)
			% We have to find Q
			% We align the nearest normal n (called normal here)
			% So we find the transforms R and P such that -
			% pT = R [1 0 0]'		T = P [1 0 0]'
			% pn = R [0 1 0]'		n = P [0 1 0]'
			% 
			% Also we find the transform Z such that
			% S-1 pn = Z [0 1 0]'
			% and S-1 pT = Z [1 0 0]' (the tangent vector doesn't change under
			%                          this transform)
			% Therefore Z = S-1 R
			%
			% Note that for alignment of normal n and tangent T, the
			% following must hold -
			% Q-1 n = Z [0 1 0]'
			% Therefore Z = Q-1 P
			% => Q-1 P = S-1 R
			% => Q = P R-1 S
			%
			R	= quatFromFrame( prevT, prevNormal );
			P	= quatFromFrame( T, normal );
			prim.q	= QuatProd( P, QuatProd( QuatInv(R), prevPrim.q ) );
		end

		% Save back.
		tube(i)		= prim;
		prevPrim	= prim;
	end
end

