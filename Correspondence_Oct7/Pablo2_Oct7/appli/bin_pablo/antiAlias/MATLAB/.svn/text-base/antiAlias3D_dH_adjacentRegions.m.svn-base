function phi = antiAlias3D_dH_adjacentRegions(img, verbose, movementThreshold, bandwidth, maxIterations, spacing)
%
% Anti-aliases by using a dH flow.
%
% function phi = antiAlias3D_K1(img)
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

import java.util.HashMap;

if( nargin < 2 )
	verbose = 2;
end

% some parameters
if( nargin < 4 )
	bandwidth		= 3;
end
if( nargin < 5 )
	maxIterations	= 5000;
end
if( nargin < 6 )
	spacing = [1 1 1];
else
	spacing = spacing ./ min(spacing);
end

if( nargin < 3 )
	% the absolute maximum movement allowed.
	% For true anti-aliasing, this should be 0.5, anything higher will
	% cause smoothing of data. Also, larger values of movementThreshold
	% should be accompanied by larger bandwidths.
	movementThreshold   = 0.5;
end


% dt			= 0.0625;		% should be less than (1/2^(#dims+1))
eps			= 0.01;			% if max(curvature) is less than this,
							% then assume that the level set
							% evolution has converged.

% Are we already passed in a distance map?
if( strcmp(class(img),'double') == 1 || strcmp(class(img),'single') == 1 )
	phi = single(img);
	if( verbose >= 1 )
		fprintf('Using passed-in distance map\n');
	end
else
	if( ~islogical(img) )
		img = logical(img);
	end
	% compute the distance map.
	phi = zeros( size( img ), 'single' );
	dist = bwdist( img, 'quasi-euclidean' ) - 0.5;	% outside
	phi( ~img ) = dist( ~img );
	dist = -bwdist( ~img, 'quasi-euclidean' ) + 0.5;	% inside
	phi( img ) = dist( img );
	% correct for the fact that we assume pixel values to be at the center of
	% the cell and the contour that is to be extracted is the one for 0.5
	% (between 0 and 1) (matlab bwdist thinks differently)
	clear dist
	if( verbose >= 1 )
		fprintf('Computed initial level set.\n');
	end
end

% Find the mask for the band of the computation.
img = -bandwidth < phi & phi < bandwidth;
% mask off all boundary regions
img([1 end], :, :) = false;
img(:, [1 end], :) = false;
img(:, :, [1 end]) = false;

% find the indices for the mask, and then find the indices for all the
% derivatives. The x and y usage here is consistent with ndgrid instead of
% meshgrid, but that doesn't matter.
i000 = uint32(find(img));

% the intial level set.
phi0 = phi(i000);

sz = size(img);
clear img;
[x,y,z] = ind2sub(sz, i000);
x = uint32(x);
y = uint32(y);
z = uint32(z);
im00 = uint32(sub2ind( sz, x-1, y  , z  ));
ip00 = uint32(sub2ind( sz, x+1, y  , z  ));
i0m0 = uint32(sub2ind( sz, x  , y-1, z  ));
i0p0 = uint32(sub2ind( sz, x  , y+1, z  ));
i00m = uint32(sub2ind( sz, x  , y  , z-1));
i00p = uint32(sub2ind( sz, x  , y  , z+1));

imm0 = uint32(sub2ind( sz, x-1, y-1, z  ));
ipp0 = uint32(sub2ind( sz, x+1, y+1, z  ));
imp0 = uint32(sub2ind( sz, x-1, y+1, z  ));
ipm0 = uint32(sub2ind( sz, x+1, y-1, z  ));

im0m = uint32(sub2ind( sz, x-1, y  , z-1));
ip0p = uint32(sub2ind( sz, x+1, y  , z+1));
im0p = uint32(sub2ind( sz, x-1, y  , z+1));
ip0m = uint32(sub2ind( sz, x+1, y  , z-1));

i0mm = uint32(sub2ind( sz, x  , y-1, z-1));
i0pp = uint32(sub2ind( sz, x  , y+1, z+1));
i0mp = uint32(sub2ind( sz, x  , y-1, z+1));
i0pm = uint32(sub2ind( sz, x  , y+1, z-1));

ordinates = [x y z];
clear x y z;

% the threshold for every vertex
upperThreshold = zeros(numel(i000), 1, 'single');
lowerThreshold = zeros(numel(i000), 1, 'single');

% build neighbor list for every vertex.
ineighbors = zeros(numel(i000)+1, 18, 'uint32');
neighCount = zeros(numel(i000)+1, 1, 'single');
% create a hash table to do a quick reverse lookup in i000.
ji000 = HashMap(numel(i000), 0.5);
for i = 1:numel(i000)
	ji000.put(i000(i), uint32(i));
end
% A temporary variable to store neighbors (make sure this is a row vecotr).
% change to size 18 if using the 18 neighbor mask.
neighs = zeros(1,18, 'uint32');
for vi = 1:numel(i000)
	% 6 or 18 neighbor mask?
%  	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi)];
 	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi) ...
		imm0(vi), ipp0(vi), imp0(vi), ipm0(vi), ...
		im0m(vi), ip0p(vi), im0p(vi), ip0m(vi), ...
		i0mm(vi), i0pp(vi), i0mp(vi), i0pm(vi) ];
	% TODO: Figure a way to do this for movementThresholds > 1
	% upperThreshold is for inside (-ve signed distance)
	% lowerThreshold is for outside (+ve signed distance)
	if( movementThreshold > 1 )
		upperThreshold(vi) = movementThreshold;
		lowerThreshold(vi) = movementThreshold;
	else
		if( any(phi(pt) >= phi0(vi) + movementThreshold ) )
			upperThreshold(vi) = movementThreshold;
		else
			upperThreshold(vi) = 0.5;
		end
		if( any(phi(pt) <= phi0(vi) - movementThreshold ) )
			lowerThreshold(vi) = movementThreshold;
		else
			lowerThreshold(vi) = 0.5;
		end
	end
	neighsCount = 0;
	for ipt = 1:numel(pt)
		% do a reverse lookup for the index of the point.
		ni = ji000.get(pt(ipt));
		if( ~isempty(ni) )
			neighsCount = neighsCount + 1;
			neighs(neighsCount) = ni;
		end
	end
	neighs(neighsCount+1:end) = numel(i000)+1;
	ineighbors(vi,:) = neighs;
	neighCount(vi) = neighsCount;
end
ineighbors(end,:) = numel(i000)+1;
clear neighs neighsCount vi pt ipt ni;

if( verbose >= 1 )
	fprintf('Computed band and its indices.\n');
end

if( verbose >= 2 )
	fig1 = figure(1);
	set(0,'CurrentFigure',fig1);
	clf;
% 	origSurface = patch(isosurface(phi, 1e-100));
% 	set(origSurface, 'FaceColor', [0 1 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
	view([331 20]);
	camlight; lighting gouraud;
	xlabel('x');
	ylabel('y');
	zlabel('z');
	axis vis3d;
	axis on;
	box on;
	daspect(1./spacing);
	p = patch(isosurface(phi,1e-100));
	set(p, 'FaceColor', 'red', 'EdgeColor', 'none');

	fig2 = figure(2);
	backdrop = zeros(size(phi),'single');

	% Added these variables here to facilitate debugging.
	drawIsoSurface = true;
	steps = 20;
	xs = ceil(size(phi,2)/2);
	ys = 1:size(phi,1);
	zs = 1:size(phi,3);
% 	xs = 420;
% 	ys = 198:320;
% 	zs = 29:53;

	initZeroLevelSet = contour(squeeze(phi(ys,xs,zs)), 1e-100);
	initZeroLevelSet = initZeroLevelSet(:,2:end);
end

% mean curvature - pre-allocated.
% the one extra element is to store a zero.
H = zeros(numel(i000)+1,1,'single');
K = H;
iNZ = false(size(H));

spacing = [1 1 1];


iter = 0;
while iter <= maxIterations
	iter = iter + 1;
	phix	= (phi(ip00) - phi(im00)) ./ ( 2*spacing(:,1) );
	phiy	= (phi(i0p0) - phi(i0m0)) ./ ( 2*spacing(:,2) );
	phiz	= (phi(i00p) - phi(i00m)) ./ ( 2*spacing(:,3) );
	phixx	= (phi(im00) - 2*phi(i000) + phi(ip00)) ./ (spacing(:,1).^2);
	phiyy	= (phi(i0m0) - 2*phi(i000) + phi(i0p0)) ./ (spacing(:,2).^2);
	phizz	= (phi(i00m) - 2*phi(i000) + phi(i00p)) ./ (spacing(:,3).^3);
	phixy	= (phi(ipp0) + phi(imm0) - phi(ipm0) - phi(imp0)) ./ (4*spacing(:,1).*spacing(:,2));
	phiyz	= (phi(i0pp) + phi(i0mm) - phi(i0pm) - phi(i0mp)) ./ (4*spacing(:,2).*spacing(:,3));
	phixz	= (phi(ip0p) + phi(im0m) - phi(ip0m) - phi(im0p)) ./ (4*spacing(:,1).*spacing(:,3));

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom	= phix_2 + phiy_2 + phiz_2;

	% non zero elements.
	iNZ(1:end-1) = denom > 0;
	% mean curvature
	H(~iNZ) = 0;
	H(iNZ) = ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
			+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
			+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
			- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
			- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
			- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

	% gaussian curvature
	K(~iNZ) = 0;
	K(iNZ) = (phix_2(iNZ).*(phiyy(iNZ).*phizz(iNZ) - phiyz(iNZ).^2) ...
			+ phiy_2(iNZ).*(phixx(iNZ).*phizz(iNZ) - phixz(iNZ).^2) ...
			+ phiz_2(iNZ).*(phixx(iNZ).*phiyy(iNZ) - phixy(iNZ).^2) ...
		+ 2*(phix(iNZ).*phiy(iNZ).*(phixz(iNZ).*phiyz(iNZ) - phixy(iNZ).*phizz(iNZ)) ...
		   + phiy(iNZ).*phiz(iNZ).*(phixy(iNZ).*phixz(iNZ) - phiyz(iNZ).*phixx(iNZ)) ...
		   + phix(iNZ).*phiz(iNZ).*(phixy(iNZ).*phiyz(iNZ) - phixz(iNZ).*phiyy(iNZ)))) ...
		./ denom(iNZ).^2;

	if( iter == 0 )
		% we will sweep along z axis.
		% this is a special iteration where we identify all the flat
		% regions and set up the anisotropic spacing for each.
		vis   = find( abs(H(1:end-1)) < 1e-5 & abs(K(1:end-1)) < 1e-5 );
		spacing = ones([numel(i000) 3]);
		x = 1;
		y = 2;
		z = 3;
		for idim = 1:3
			region  = zeros(size(vis), 'uint32');
			stack   = zeros(size(vis), 'uint32');
			visited = false(size(i000));
			flatRegions = containers.Map(double(0), uint32([0 0]));
			remove(flatRegions, 0);
			regions = cell(size(vis));
			numRegions = 0;
			for vi = vis'
				if( visited(vi) )
					continue;
				end
				level = floor(0.5+phi0(vi));
				zi = ordinates(vi,z);
				npts        = 1;
				stack(npts) = vi;
				visited(vi) = true;
				regionSize  = 0;
				while( npts > 0 )
					pt   = stack(npts);
					npts = npts - 1;
					regionSize = regionSize + 1;
					region(regionSize) = pt;
					% look at neighbors and add them if they are part of the
					% region.
					neighs = ineighbors(pt,:);
					neighs(neighs == numel(i000) + 1) = [];
					% part of connected flat region, push the point(s) on
					% the stack.
					neighs = neighs(ordinates(neighs, z) == zi ...
							& floor(0.5+phi0(neighs)) == level & ~visited(neighs) ...
							& abs(H(neighs)) < 1e-5 & abs(K(neighs)) < 1e-5 );
					visited(neighs)   = true;
					npts = npts + numel(neighs);
					stack(npts-numel(neighs)+1:npts) = neighs;
				end
				% Include all flat regions. Singly thick regions will be
				% eliminated later through opening operations.
				numRegions = numRegions + 1;
				regions{numRegions} = region(1:regionSize)';
				if( isKey(flatRegions, zi) )
					flatRegions(zi) = [flatRegions(zi) numRegions];
				else
					flatRegions(zi) = numRegions;
				end
			end
			base = zeros(sz(x), sz(y));
			% now that we have list of all regions, sweep through the list
			for zi = 1:sz(z)
				% get all points in the flat regions in this plane.
				if( ~isKey(flatRegions, zi) )
					% just ignore and go on if there are none.
					continue;
				end
				if( verbose >= 1 )
					fprintf('Processing ordinate %d\n', zi);
				end
				if(isKey(flatRegions,zi-1) && isKey(flatRegions,zi+1))
					adjacentRegions = [flatRegions(zi-1) flatRegions(zi+1)];
				elseif( isKey(flatRegions,zi-1) )
					adjacentRegions = flatRegions(zi-1);
				elseif( isKey(flatRegions,zi+1) )
					adjacentRegions = flatRegions(zi+1);
				else
					% no adjacent flat regions? hmm, no anisotropic spacing
					% here
					continue;
				end
				% now process for each level individually.
				for level = -bandwidth:bandwidth
					levAdjacent = true(size(adjacentRegions));
					for ri = 1:numel(adjacentRegions)
						if( floor(0.5+phi0(regions{adjacentRegions(ri)}(1))) ~= level )
							levAdjacent(ri) = false;
						end
					end
					levAdjacent = adjacentRegions(levAdjacent);
					if( isempty(levAdjacent) )
						continue;
					end
					% The following checks for containment
					% of the below and upper regions and process regions in this
					% ordinate one at a time, ignoring below and upper regions that
					% are outside the region.
					% Generate a binary image with adjacent flat regions.
					slice = base;
					slice(sub2ind(sz([x,y]), ordinates([regions{levAdjacent}], x), ordinates([regions{levAdjacent}], y) )) = 1;
					% fill any holes in this image.
					adjacent = imfill(slice, 'holes');
					% now loop over all current regions.
					for region = flatRegions(zi)
						if( floor(0.5+phi0(regions{region}(1))) ~= level )
							% process for the same levels. (Do not mix one
							% level with another.)
							continue;
						end
						slice = base;
						pts = regions{region};
						slice(sub2ind(sz([x,y]), ordinates(pts, x), ordinates(pts, y) )) = 1;
						slice = imfill(slice, 'holes');
						slice = imopen(slice, strel(1));
						slice = slice & adjacent;
						if( ~any(slice(:)) )
							% if there is no intersection with this region,
							% then there's nothing to process. TODO: we should
							% set a spacing relative to some measure of area.
							% This could be computed via some medial transform.
							continue;
						end
						% Compute distance map of this image. Note, these distances
						% should be from pixel centers. Therefore matlab's bwdist works
						% correctly in this case.
						[slice, closestPoint] = bwdist(slice);
						% Pick up spacing for each point and move along lines
						% to the closest boundary point setting spacing for all
						% other points accordingly.
						for pt = double([ordinates(pts,x)'; ordinates(pts,y)'])
							% find the closest boundary point.
							[tox, toy] = ind2sub(sz([x,y]), closestPoint(pt(1), pt(2)) );
							to = [tox; toy];
							length = sqrt(sum((pt-to).^2));
							for t = linspace(0,1,length+1);
								tmp = floor(pt.*t + to.*(1-t));
								ipt = [0 0 0];
								ipt(x) = tmp(1);
								ipt(y) = tmp(2);
								ipt(z) = zi;
								ipt = ji000.get( uint32(sub2ind(sz, ipt(1), ipt(2), ipt(3))) );
								if( ~isempty(ipt) )
									spacing(ipt, z) = min( spacing(ipt,z), 1/slice(pt(1), pt(2)) );
								end
							end
						end
					end
				end
			end
			tmp = circshift([x;y;z], 1);
			x = tmp(1);
			y = tmp(2);
			z = tmp(3);
		end
	end

	discriminant = H.^2 - 2*K;
	sel = discriminant > 0;
	K1 = zeros(size(H), 'single');
	K2 = K1;
	K1(sel) = (H(sel) + sqrt(discriminant(sel)))*0.5;
	K2(sel) = (H(sel) - sqrt(discriminant(sel)))*0.5;
% 	C = sqrt((K1.^2 + K2.^2)/2);
% 	%C = 2*log(C)/pi
% 	S = K1;
% 	sel = (K1~=K2);
% 	S(sel) = 2*atan( (K1(sel)+K2(sel))./(K2(sel)-K1(sel)) ) / pi;
% 	S(~sel) = sign(K1(~sel));

	dH = H - sum(H(ineighbors),2)./neighCount;
	dK1 = K1 - sum(K1(ineighbors),2)./neighCount;
	dK2 = K2 - sum(K2(ineighbors),2)./neighCount;
% 	dS = S - sum(S(ineighbors),2)./neighCount;

	% convex points do not flow out and concave points do not flow in --
	% sign conventions of distance map is important for this to work
	% correctly!
	update = zeros(numel(H)-1,1,'single');
	sel = K(1:end-1) < 0;
	update(sel) = dH([sel; false]);
	sel = K(1:end-1) >= 0 & H(1:end-1) > 0;
	update(sel) = max(dK1([sel; false]),0);
% 	update(sel) = dK1([sel; false]);
	sel = K(1:end-1) >= 0 & H(1:end-1) < 0;
	update(sel) = min(dK2([sel; false]),0);
% 	update(sel) = dK2([sel; false]);

	% restricted mean curvature flow.
% 	update = zeros(numel(H)-1,1,'single');
% 	sel = H(1:end-1) >= 0;
% 	update(sel) = max(dK1([sel; false]),0);
% 	sel = H(1:end-1) < 0;
% 	update(sel) = min(dK2([sel; false]),0);

	% S weighted dK1-dK2 flow.
% 	update = (dK2.*(S+1) + dK1.*(S-1)) / 2;
% 	update(end) = [];

% 	% relative movement weighted dK1-dK2 flow.
% 	update = (dK1.*abs(dK1./K1) + dK2.*abs(dK2./K2))./(abs(dK1./K1) + abs(dK2./K2));
% 	update(end) = [];
	
	update = update.*(denom.^0.5);
	
	% restrict update to a maximum of some value so that the time step
	% doesn't become ridiculously small.
	update = max(min(update, 10), -10);

	error = max(abs(update));
	if( error < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	% Compute dt adaptively with each iteration.
	% The 0.01 in the denominator is to ensure that dt doesn't grow to
	% infinity.
	dt = 1/max(max(abs(update)),16);
	if( verbose >= 1 )
		fprintf('\n%4d: %6.4f, %6.4f, %f', iter, error, sqrt(mean(update.^2)), dt );
	end

	phi(i000) = max(min(phi(i000) + dt*update, phi0 + upperThreshold), phi0 - lowerThreshold);

	if( mod(iter, steps) == 0 )
		dumpStatus;
	end
end
if( verbose >= 1 )
	fprintf('\n');
end

fprintf('Last maximum error %d\n', max(abs(update)) );

	function dumpStatus
		if( verbose >= 2 )
			if( drawIsoSurface )
				set(0,'CurrentFigure',fig1);
				delete(p);
				p = patch(isosurface(phi,1e-100));
				set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
			end

			set(0,'CurrentFigure',fig2);
			clf;

			subplot(2,2,1);
			backdrop(i000) = H(1:end-1);
			imagesc(squeeze(backdrop(ys,xs,zs)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			c = contour(squeeze(phi(ys,xs,zs)), 1e-100,'k');
			c = c(:,2:end);

			subplot(2,2,2);
			backdrop(i000) = 1./spacing(:,2);
			imagesc(squeeze(backdrop(ys,xs,zs)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,3);
			backdrop(i000) = K(1:end-1);
			imagesc(squeeze(backdrop(ys,xs,zs)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,4);
			backdrop(i000) = update;
			imagesc(squeeze(backdrop(ys,xs,zs)), [-1 1]);
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			drawnow;
		end
	end
end
