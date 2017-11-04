function [syms]	= minus( aprims, bprims )
%
% sym = aprims - bprims
% 
% Overloads the subtraction operator for tube sym primitives.
% The result is in manifold space and a and b should also be in manifold
% space.
%

if( sum([aprims.inTangentSpace bprims.inTangentSpace]) > 0 )
	disp('aprims = ' );
	aprims
	disp('bprims = ' );
	bprims
	error('Parameters aprims and bprims should be in manifold space.');
end

syms	= aprims;

for i = 1:prod(size(syms))
	a	= aprims(i);
	b	= bprims(i);
	sym	= a;
	sym.pos			= a.pos - b.pos;
	sym.r			= a.r/b.r;
% TODO: Temporary FIX
%	sym.dr			= a.dr - b.dr;
	sym.dr			= a.dr;
	% TODO: verify RP1 stats.
	sym.hca			= atan(tan(a.hca-pi/2.0) - tan(b.hca-pi/2.0)) + pi/2.0;
	sym.elongation	= a.elongation / b.elongation;
	sym.U0			= symAction( b.U0, a.U0 );
	sym.inTangentSpace	= 0;
	syms(i)	= sym;
end

return;
