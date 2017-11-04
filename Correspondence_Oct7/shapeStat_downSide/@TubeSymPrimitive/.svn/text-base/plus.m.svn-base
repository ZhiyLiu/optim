function [tprims] = plus( aprims, bprims );
%
% [t] = a + b;
% 
% Adds the two symmetric atoms a and b
% (composition)
%
% The result is in manifold space if a and b are in manifold
% space and vice-versa
%
% a and b can be arrays of primitives.
%


if( length(aprims) ~= length(bprims) )
	error('Parameters a and b should be of the same size');
end

tprims	= aprims;

if( (sum([aprims.inTangentSpace]) == length(aprims)) && (sum([bprims.inTangentSpace]) == length(bprims)) )
	for i = 1:prod(size(aprims))
		a	= aprims(i);
		b	= bprims(i);
		t	= a;

		t.pos	= a.pos + b.pos;
		t.r		= a.r + b.r;
% TODO: Temporary FIX
%		t.dr	= a.dr + b.dr;
		t.dr	= a.dr;
		t.elongation	= a.elongation + b.elongation;
		% TODO: verify RP1 stats.
		t.hca	= a.hca + b.hca;
		t.U0	= SphereLog(QuatRotVec( QuatInv(getRotation(SphereExp(a.U0))), SphereExp(b.U0) ));

		tprims(i)	= t;
	end
elseif( (sum([aprims.inTangentSpace]) == 0) && (sum([bprims.inTangentSpace]) == 0) )
	for i = 1:prod(size(aprims))
		a	= aprims(i);
		b	= bprims(i);
		t	= a;

		t.pos	= a.pos + b.pos;
		t.r		= a.r * b.r;
% TODO: Temporary FIX
%		t.dr	= a.dr + b.dr;
		td.r	= a.dr;
		t.elongation	= a.elongation .* b.elongation;
		% TODO: verify RP1 stats.
		t.hca	= atan( tan(a.hca-pi/2.0) + tan(b.hca-pi/2.0) ) + pi/2.0;
		t.U0	= QuatRotVec( QuatInv(getRotation(a.U0)), b.U0 );

		tprims(i)	= t;
	end
else
	disp('a = ' );
	aprims
	disp('b = ' );
	bprims
	error('Parameters a and b should be both in manifold space or both in tangent space.');
end

return;
