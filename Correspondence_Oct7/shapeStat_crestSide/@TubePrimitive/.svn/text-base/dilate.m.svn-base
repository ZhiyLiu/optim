function [prims] = dilate(prims, amount )
%
% function [prims] = dilate(prims, amount)
%
% Dilates the primitives by the specified amount.
% Dilation is simply addition of amount to the radius.
%

[r, c] = size(prims);

for i=1:r
	for j=1:c
		p	= prims(i,j);
		p.r	= p.r + amount;
		prims(i,j)	= p;
	end
end
