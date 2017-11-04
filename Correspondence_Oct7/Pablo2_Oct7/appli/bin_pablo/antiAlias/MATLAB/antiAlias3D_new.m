function phi = antiAlias3D(img, verbose)
%
% function phi = antiAlias3D(img)
%	Anti aliases a 3D binary image (1 inside, 0 outside and the boundary
% assumed to lie at 0.5); return a level-set with 0 representing the
% boundary.
%
% Input:
%	img		The input binary image (a 3D array)
%	verbose	0 for no progress, 1 for text, 2 for text+visuals (default).
% Output:
%	phi		The anti-aliased level set
%
% Notes:
%	Increasing the bandwidth or changing the method for computing the
% distance map from euclidean to quasi-euclidean doesn't significantly
% affect the results.
%

import java.util.LinkedList;
import java.util.HashMap;

if( nargin < 2 )
	verbose = 2;
end

% some parameters
maxIterations		= 1000;
bandwidth			= 4;
constraintLeeway	= 0.3;	% by how much should we allow the hard
							% constraints to be disrespected.

% .. and other not so parameter-like constants.
dt					= 0.0625;		% should be less than (1/2^(#dims+1))
movementThreshold	= sqrt(3)/2;	% the absolute maximum movement
									% allowed. This is sqrt(#dims)/2 i.e.
									% length of body diagonal of a voxel.
eps = 1e-2;							% if max(curvature) is less than this,
									% then assume that the level set
									% evolution has converged.

% compute the distance map.
% dout = bwdist( img );
% din = -bwdist( 1-img );
dout = bwdist( img, 'quasi-euclidean' );
din = -bwdist( 1-img, 'quasi-euclidean' );
% correct for the fact that we assume pixel values to be at the center of
% the cell and the contour that is to be extracted is the one for 0.5
% (between 0 and 1) (matlab bwdist thinks differently)
dout = dout - 0.5;
din  = din + 0.5;
% put them together to one distance map
mask = logical(img);
phi = zeros( size( mask ), 'single' );
phi( mask ) = din( mask );
phi( ~mask ) = dout( ~mask );
clear dout din
if( verbose >= 1 )
	fprintf('Computed initial level set.\n');
end

% Find the mask for the band of the computation.
mask = -bandwidth < phi & phi < bandwidth;
% mask off all boundary regions
mask([1 end], :, :) = false;
mask(:, [1 end], :) = false;
mask(:, :, [1 end]) = false;

% We will use a separate array to disable computations for those regions
% that have reached the stopping criteria for evolution to avoid
% recomputing the indices array -- evolveMask.

% find the indices for the mask, and then find the indices for all the
% derivatives. The x and y usage here is consistent with ndgrid instead of
% meshgrid, but that doesn't matter.
i000 = find(mask);

% create a hash table to do a quick reverse lookup in i000.
ji000 = HashMap(numel(i000));
for i = 1:numel(i000)
	ji000.put(i000(i), i);
end

sz = size(mask);
clear mask;
[x,y,z] = ind2sub(sz, i000);
im00 = sub2ind( sz, x-1, y  , z  );
ip00 = sub2ind( sz, x+1, y  , z  );
i0m0 = sub2ind( sz, x  , y-1, z  );
i0p0 = sub2ind( sz, x  , y+1, z  );
i00m = sub2ind( sz, x  , y  , z-1);
i00p = sub2ind( sz, x  , y  , z+1);

imm0 = sub2ind( sz, x-1, y-1, z  );
ipp0 = sub2ind( sz, x+1, y+1, z  );
imp0 = sub2ind( sz, x-1, y+1, z  );
ipm0 = sub2ind( sz, x+1, y-1, z  );

im0m = sub2ind( sz, x-1, y  , z-1);
ip0p = sub2ind( sz, x+1, y  , z+1);
im0p = sub2ind( sz, x-1, y  , z+1);
ip0m = sub2ind( sz, x+1, y  , z-1);

i0mm = sub2ind( sz, x  , y-1, z-1);
i0pp = sub2ind( sz, x  , y+1, z+1);
i0mp = sub2ind( sz, x  , y-1, z+1);
i0pm = sub2ind( sz, x  , y+1, z-1);

clear x y z sz;

if( verbose >= 1 )
	fprintf('Computed band and its indices.\n');
end

if( verbose >= 2 )
	figure
	view([-37.5 30]);
end

% the intial level set.
phi0 = phi(i000);

for flowDir = [-1 1]
% shrink image by half a voxel (this is to account for the fact that we
% only allow our level-set to grow and not shrink.)
% Do not do this if you want to try out the min-max flow.
% 0.5 is half of movementThreshold.
% if( flowDir == -1 )
% 	% erode by half a voxel.
% 	phi  = phi + 0.5;
% 	phi0 = phi0 + 0.5;
% else
% 	% dilate by half a voxel (this will dilate existing expanded regions
% 	% from previous min flow too :( )
% 	phi  = phi - 0.5;
% 	phi0 = phi0 - 0.5;
% end

evolveMask = true(size(i000));
adaptiveDt = repmat(dt, [numel(phi0), 1]);
if( verbose >= 1 )
	fprintf('Iterations     ');
end
for iter = 0:maxIterations
	if( verbose >= 1 )
		fprintf('\b\b\b\b%4d', iter);
	end
	if( mod(iter, (maxIterations/100)) == 0 && verbose >= 2)
		[az, el] = view;
		clf;
		isosurface(phi,0);
		axis vis3d;
		axis on;
		box on;
		axis equal;
		view(az,el);
		drawnow;
		pause(0.1);
	end
	phix	= (phi(ip00) - phi(im00)) * 0.5;
	phiy	= (phi(i0p0) - phi(i0m0)) * 0.5;
	phiz	= (phi(i00p) - phi(i00m)) * 0.5;
	phixx	= phi(im00) - 2*phi(i000) + phi(ip00);
	phiyy	= phi(i0m0) - 2*phi(i000) + phi(i0p0);
	phizz	= phi(i00m) - 2*phi(i000) + phi(i00p);
	phixy	= (phi(ipp0) + phi(imm0) - phi(ipm0) - phi(imp0)) * 0.25;
	phiyz	= (phi(i0pp) + phi(i0mm) - phi(i0pm) - phi(i0mp)) * 0.25;
	phixz	= (phi(ip0p) + phi(im0m) - phi(ip0m) - phi(im0p)) * 0.25;

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom	= phix_2 + phiy_2 + phiz_2;

	H = zeros(size(denom));
	K = H;
	% non zero elements.
	iNZ = denom > 0;
	% mean curvature
	H(iNZ) = ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
			+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
			+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
			- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
			- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
			- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

	% gaussian curvature
	K(iNZ) = (phix_2(iNZ).*(phiyy(iNZ).*phizz(iNZ) - phiyz(iNZ).^2) ...
			+ phiy_2(iNZ).*(phixx(iNZ).*phizz(iNZ) - phixz(iNZ).^2) ...
			+ phiz_2(iNZ).*(phixx(iNZ).*phiyy(iNZ) - phixy(iNZ).^2) ...
		+ 2*(phix(iNZ).*phiy(iNZ).*(phixz(iNZ).*phiyz(iNZ) - phixy(iNZ).*phizz(iNZ)) ...
		   + phiy(iNZ).*phiz(iNZ).*(phixy(iNZ).*phixz(iNZ) - phiyz(iNZ).*phixx(iNZ)) ...
		   + phix(iNZ).*phiz(iNZ).*(phixy(iNZ).*phiyz(iNZ) - phixz(iNZ).*phiyy(iNZ)))) ...
		./ denom(iNZ).^2;
	% normal
	n = zeros(numel(denom),3);
	n(iNZ,:) = [ phix(iNZ), phiy(iNZ), phiz(iNZ) ] ./ repmat(denom(iNZ).^0.5, [1 3]);
	n(iNZ,:) = 0.5 * n(iNZ,:) ./ repmat( max(abs(n(iNZ,: )),[],2), [1 3] );

	movementConstraint = 0.5*ones(size(denom));
	movementConstraint(iNZ) = sqrt(sum(n(iNZ,:).^2,2));
	movementConstraint = min(movementConstraint + constraintLeeway,movementThreshold);
	% the min-max flow -- bad.
%	C = zeros(size(H));
% 	C(img(i000) == 1) = min(H(img(i000) == 1),0);
% 	C(img(i000) == 0) = max(H(img(i000) == 0),0);

    if( all(K >= 0 & H >= 0) )
        % This is an extremely rare condition to reach. It prevents an
        % ellipsoid from collapsing.
		fprintf('\nLevel set has been convexified.');
		break;
    end

	% the min-flow.
	if( flowDir < 0 )
		C = min(H, 0);
	else
		C = max(H, 0);
	end

	% flow all concave (H < 0, K > 0), hyperbolic (K < 0),
	% and parabolic interface between concave and hyperbolic regions
	% (H < 0, K == 0) using mean curvature.
% 	C = zeros(size(H));
% 	sel	   = (K > 0 & H < 0) | (K < 0) | (K == 0 & H < 0);
% 	C(sel) = H(sel);


	% flow only the convex regions -- sucks.
% 	C = zeros(size(H));
% 	sel	   = (K > 0 & H < 0);
% 	C(sel) = H(sel);

	% the K2 min-flow
% 	discriminant = H.^2 - 2*K;
% 	sel = discriminant > 0;
% 	K2 = zeros(size(H));
% 	K2(sel) = (H(sel) - sqrt(discriminant(sel)))*0.5;
% 	C = min(K2,0);

	if( sum(evolveMask) == 0 || max(abs(C(evolveMask))) < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end

	newphi = phi(i000) + adaptiveDt.*C;
	% Any pixel that has been displaced by more than a voxel (threshold)
	% should stop evolving. Any connecting region of this pixel,
	% connectivity given by K <= 0 will also stop evolving.
% 	elements = find( (abs(newphi - phi0) > movementThreshold) & evolveMask)';
 	elements = find( (abs(newphi - phi0) > movementConstraint) & evolveMask)';
	visited  = false( size(evolveMask) );
	for ei = elements;
		% if this pixel has not already been masked out by another point in
		% this loop, do a flood fill.
		if( evolveMask(ei) && ~visited(ei) )
			% find all connected regions with curvature less than 0 and set
			% the mask to false.
			pts = java.util.LinkedList();
			pts.addLast(ei);
			region = [];
			dtRatio = 1;
			visited(ei)   = true;
			while( ~(pts.isEmpty()) )
				pt = pts.removeFirst();
				region(end+1) = pt;
				newDt = (movementConstraint(pt) - abs(phi0(pt) - phi(i000(pt)))) / abs(C(pt));
				if( newDt < eps )
					dtRatio = 0;
				else
					% the divide by 2 is to make it the flow a little
					% conservative so that we don't overshoot.
					dtRatio = min(newDt/adaptiveDt(pt)/2, dtRatio);
				end
				% look at neighbors and add them if they are part of the
				% region. (6 neighborhood)
				pt = [im00(pt), ip00(pt), i0m0(pt), i0p0(pt), i00m(pt), i00p(pt)];
				for ipt = 1:6
					% do a reverse lookup for the index of the point.
					i = ji000.get(pt(ipt));
					if( ~isempty(i) )
						% i is now index in i000 array for the pt.
						if( ~visited(i) && evolveMask(i) && abs(phi0(i) - newphi(ei)) <= 2 ...
							 && C(i)*flowDir >= 0 )
							visited(i)   = true;
							% part of connected region, push the point in
							% the queue.
							pts.addLast(i);
						end
					end
				end
			end
			if( dtRatio == 0 )
				evolveMask(region) = false;
			else
				adaptiveDt(region) = dtRatio * adaptiveDt(region);
			end
		end
	end
	phi(i000(evolveMask)) = newphi(evolveMask);
end
if( verbose >= 1 )
	fprintf('\n');
end
end

fprintf('Last maximum error %d\n', max(abs(C(evolveMask))) );
