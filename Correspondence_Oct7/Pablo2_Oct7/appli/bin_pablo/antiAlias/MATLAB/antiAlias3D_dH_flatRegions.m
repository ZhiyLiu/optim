function [phi, Kfull, Hfull] = antiAlias3D_dH_flatRegions(img, verbose, movementThreshold, bandwidth, maxIterations, spacing)
%
% Anti-aliases by using a dH flow.
%
% function phi = antiAlias3D_dH(img)
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

dt  = 0.0625;
eps	= 0.01;			% if max(curvature) is less than this,
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

immm = uint32(sub2ind( sz, x-1, y-1, z-1));
immp = uint32(sub2ind( sz, x-1, y-1, z+1));
impm = uint32(sub2ind( sz, x-1, y+1, z-1));
impp = uint32(sub2ind( sz, x-1, y+1, z+1));
ipmm = uint32(sub2ind( sz, x+1, y-1, z-1));
ipmp = uint32(sub2ind( sz, x+1, y-1, z+1));
ippm = uint32(sub2ind( sz, x+1, y+1, z-1));
ippp = uint32(sub2ind( sz, x+1, y+1, z+1));

allNeighbors = [ im00 ip00 i0m0 i0p0 i00m i00p ...
	imm0 ipp0 imp0 ipm0 ...
	im0m ip0p im0p ip0m ...
	i0mm i0pp i0mp i0pm ...
	immm immp impm impp ipmm ipmp ippm ippp ];

ordinates = [x y z];
clear x y z;

% the threshold for every vertex
upperThreshold = zeros(numel(i000), 1, 'single');
lowerThreshold = zeros(numel(i000), 1, 'single');

% build neighbor list for every vertex.
ineighbors = zeros(numel(i000)+1, 6, 'uint32');
% direction of the nieghboring element: 1 along x, 2 along y, 3 along z.
idirs      = zeros(numel(i000)+1, 6, 'uint32');
neighCount = zeros(numel(i000)+1, 1, 'single');
% create a hash table to do a quick reverse lookup in i000.
ji000 = HashMap(numel(i000), 0.5);
for i = 1:numel(i000)
	ji000.put(i000(i), uint32(i));
end
% A temporary variable to store neighbors (make sure this is a row vecotr).
% change to size 18 if using the 18 neighbor mask.
neighs = zeros(1,6, 'uint32');
dirs   = zeros(1,6, 'uint32');
dir    = [1, 1, 2, 2, 3, 3];
for vi = 1:numel(i000)
	% 6 or 18 neighbor mask?
 	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi)];
%  	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi) ...
% 		imm0(vi), ipp0(vi), imp0(vi), ipm0(vi), ...
% 		im0m(vi), ip0p(vi), im0p(vi), ip0m(vi), ...
% 		i0mm(vi), i0pp(vi), i0mp(vi), i0pm(vi) ];
	% TODO: Figure a way to do this for movementThresholds > 1
	% upperThreshold is for inside (-ve signed distance)
	% lowerThreshold is for outside (+ve signed distance)
% 	if( movementThreshold > 1 )
		upperThreshold(vi) = movementThreshold(1);	% inside
		lowerThreshold(vi) = movementThreshold(2);	% outside
% 	else
% 		if( any(phi(allNeighbors(vi,:)) >= phi0(vi) + movementThreshold ) )
% 			upperThreshold(vi) = movementThreshold;
% 		else
% 			upperThreshold(vi) = 0;
% 		end
% 		if( any(phi(allNeighbors(vi,:)) <= phi0(vi) - movementThreshold ) )
% 			lowerThreshold(vi) = movementThreshold;
% 		else
% 			lowerThreshold(vi) = 0;
% 		end
% 	end
	neighsCount = 0;
	for ipt = 1:numel(pt)
		% do a reverse lookup for the index of the point.
		ni = ji000.get(pt(ipt));
		if( ~isempty(ni) )
			neighsCount = neighsCount + 1;
			neighs(neighsCount) = ni;
			dirs(neighsCount) = dir(ipt);
		end
	end
	neighs(neighsCount+1:end) = numel(i000)+1;
	dirs(neighsCount+1:end)   = 1;
	ineighbors(vi,:) = neighs;
	idirs(vi,:)      = dirs;
	neighCount(vi) = neighsCount;
end
ineighbors(end,:) = numel(i000)+1;
idirs(end,:) = 1;
clear neighs neighsCount vi pt ipt ni;

if( verbose >= 1 )
	fprintf('Computed band and its indices.\n');
end

if( verbose >= 2 )
	fig1 = figure(1);
	set(0,'CurrentFigure',fig1);
	clf;
	view([331 20]);
	xlabel('x');
	ylabel('y');
	zlabel('z');
	axis vis3d;
	axis on;
	box on;
	colormap(hot);
	colorbar;
	daspect(1./spacing);
% 	origSurface = patch(isosurface(phi, 1e-100));
% 	set(origSurface, 'FaceColor', [0 1 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
 	camlight; lighting gouraud;

	fig2 = figure(2);
	backdrop  = zeros(size(phi),'single');

	% Added these variables here to facilitate debugging.
	drawIsoSurface = true;
 	xs = ceil(size(phi,2)/2);
	ys = 1:size(phi,1);
 	zs = 1:size(phi,3);

% 	xs = 420;
% 	ys = 198:320;
% 	zs = 29:53;
% 
% 	xs = 299;
% 	ys = 300:500;
% 	zs = 70:90;

	initZeroLevelSet = contour(squeeze(phi(ys,xs,zs)), 1e-100);
	initZeroLevelSet = initZeroLevelSet(:,2:end);
end
steps = 21;

% mean curvature - pre-allocated.
% the one extra element is to store a zero.
H = zeros(numel(i000)+1,1,'single');
K = H;
iNZ = false(size(H));

spacing = ones([numel(H) 3]);
distances = ones(size(ineighbors));
update = zeros(size(i000));

for iter = 1:maxIterations
	phix	= (phi(ip00) - phi(im00)) ./ ( 2*spacing(1:end-1,1) );
	phiy	= (phi(i0p0) - phi(i0m0)) ./ ( 2*spacing(1:end-1,2) );
	phiz	= (phi(i00p) - phi(i00m)) ./ ( 2*spacing(1:end-1,3) );
	phixx	= (phi(im00) - 2*phi(i000) + phi(ip00)) ./ (spacing(1:end-1,1).^2);
	phiyy	= (phi(i0m0) - 2*phi(i000) + phi(i0p0)) ./ (spacing(1:end-1,2).^2);
	phizz	= (phi(i00m) - 2*phi(i000) + phi(i00p)) ./ (spacing(1:end-1,3).^3);
	phixy	= (phi(ipp0) + phi(imm0) - phi(ipm0) - phi(imp0)) ./ (4*spacing(1:end-1,1).*spacing(1:end-1,2));
	phiyz	= (phi(i0pp) + phi(i0mm) - phi(i0pm) - phi(i0mp)) ./ (4*spacing(1:end-1,2).*spacing(1:end-1,3));
	phixz	= (phi(ip0p) + phi(im0m) - phi(ip0m) - phi(im0p)) ./ (4*spacing(1:end-1,1).*spacing(1:end-1,3));

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom	= phix_2 + phiy_2 + phiz_2;

	% non zero elements.
	iNZ(1:end-1) = denom > 0.2;
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
		% some initializations for further runs.
		oldphi = phi(i000);

		% we will sweep along z axis.
		% this is a special iteration where we identify all the flat
		% regions and set up the anisotropic spacing for each.
		vis   = find( abs(H(1:end-1)) < 1e-5 & abs(K(1:end-1)) < 1e-5 );
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
					% There's always a single pixel wide hyperbolic region
					% surrounding every flat region. We will include those
					% pixels too in this flat region, but we won't add them
					% to the search stack.
					neighs = ineighbors(pt,:);
					neighs(neighs == numel(i000) + 1) = [];
					neighs = neighs(ordinates(neighs, z) == zi ...
							& floor(0.5+phi0(neighs)) == level & ~visited(neighs) );
					visited(neighs)   = true;
					regionSize = regionSize + numel(neighs);
					region(regionSize-numel(neighs)+1:regionSize) = neighs;
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
			base = false(sz(x), sz(y));
			% now that we have list of all regions, sweep through the list
			for zi = 1:sz(z)
				if( ~isKey(flatRegions, zi) )
					% just ignore and go on if there are no regions at this
					% zindex.
					continue;
				end
				if( verbose >= 1 )
					fprintf('Processing ordinate %d:%d\n', z, zi);
				end
				% get all points in the flat regions in this plane that
				% belong to the same level.
				for level = -bandwidth:bandwidth
					includeRegion = false(size(flatRegions(zi)));
					% loop over all current regions finding regions of same
					% level.
					thisRegions = flatRegions(zi);
					for ri = 1:numel(thisRegions)
						if( floor(0.5+phi0(regions{thisRegions(ri)}(1))) == level )
							% process for the same levels. (Do not mix one
							% level with another.)
							includeRegion(ri) = true;
						end
					end
					% get all points for these regions.
					pts = [regions{thisRegions(includeRegion)}];
					if( isempty(pts) )
						% if we found no regions for this level, then go
						% onto the next level.
						continue;
					end
					plane = base;
					% build an image from these regions.
					plane(sub2ind(sz([x,y]), ordinates(pts, x), ordinates(pts, y) )) = true;
					% Compute internal distance map of this image. Note, these distances
					% should be from pixel centers. Therefore matlab's bwdist works
					% correctly in this case.
					for invert = [false true]
						% This may result in asymmetry in regions which have
						% odd thickness. The center line is close towards both
						% the sides. So we should flip the plane, compute
						% the closestPoint transform and redraw the lines for
						% points which the closestPoint changes.
						if( invert )
							[plane, closestPoint] = bwdist(~plane(end:-1:1,end:-1:1));
							plane = plane(end:-1:1, end:-1:1);
							[closestPointx, closestPointy] = ind2sub( sz([x y]), closestPoint );
							closestPointx = sz(x) - closestPointx + 1;
							closestPointy = sz(y) - closestPointy + 1;
							closestPoint = sub2ind( sz([x y]), closestPointx, closestPointy );
							closestPoint = closestPoint(end:-1:1, end:-1:1);
						else
							[plane, closestPoint] = bwdist(~plane);
						end
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
								if( ~plane(tmp(1), tmp(2)) )
									continue;
								end
								ipt = [0 0 0];
								ipt(x) = tmp(1);
								ipt(y) = tmp(2);
								ipt(z) = zi;
								ipt = ji000.get( uint32(sub2ind(sz, ipt(1), ipt(2), ipt(3))) );
								if( ~isempty(ipt) )
									assert(plane(pt(1),pt(2)) >= 1);
									spacing(ipt,z) = max( spacing(ipt,z), plane(pt(1), pt(2)) );
								end
							end
						end
					end
					% now displace all points by ratio of
					% 0.5*(spacing-distance)/spacing
					% (The 0.5 is replaced by half of phi difference at
					% boundaries.)
					for ipt = pts;
						pt = double([ordinates(ipt,x); ordinates(ipt,y)]);
						from = [0 0 0];
						from(x) = pt(1);
						from(y) = pt(2);
						from(z) = zi;
						ifrom = ji000.get( uint32(sub2ind(sz, from(1), from(2), from(3))) );
						[tox, toy] = ind2sub(sz([x,y]), closestPoint(pt(1), pt(2)) );
						to = [0 0 0];
						to(x) = tox;
						to(y) = toy;
						to(z) = zi;
						phi(i000(ipt)) = phi(i000(ipt)) + ...
							(spacing(ifrom,z) - plane(pt(1), pt(2)))/spacing(ifrom,z) * ...
							(phi(to(1),to(2),to(3)) - phi(from(1),from(2),from(3)))*0.5;
					end
				end
			end
			tmp = circshift([x;y;z], 1);
			x = tmp(1);
			y = tmp(2);
			z = tmp(3);
		end
		distances = reshape(spacing( sub2ind(size(spacing), ineighbors(:), idirs(:) ) ), size(ineighbors) );
		for vi = 1:size(distances,1)
			neighCount(vi) = 0;
			for ni = 1:6
				if( ineighbors(vi,ni) ~= numel(i000) + 1 )
					neighCount(vi) = neighCount(vi) + 1/spacing(ineighbors(vi,ni), idirs(vi,ni));
				end
			end
		end
		dumpStatus;
		continue;
	end
	
	discriminant = H.^2 - 2*K;
	sel = discriminant >= 0;
	assert( all(sel) );
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

	dH = H - sum(H(ineighbors)./distances,2)./neighCount;
	dK1 = K1 - sum(K1(ineighbors)./distances,2)./neighCount;
	dK2 = K2 - sum(K2(ineighbors)./distances,2)./neighCount;
% 	dS = S - sum(S(ineighbors),2)./neighCount;

	% convex points do not flow out and concave points do not flow in --
	% sign convention of distance map is important for this to work
	% correctly! If we enable this condition, then it can create biases for
	% convex regions, in that, we would never let a convex region grow
	% outwards.
	update = zeros(numel(H)-1,1,'single');
	sel = K(1:end-1) < 0;
	update(sel) = dH([sel; false]);
	sel = K(1:end-1) >= 0 & H(1:end-1) > 0;
	update(sel) = max(dK1([sel; false]),0);
% 	update(sel) = dK1([sel; false]);
	sel = K(1:end-1) >= 0 & H(1:end-1) < 0;
	update(sel) = min(dK2([sel; false]),0);
% 	update(sel) = dK2([sel; false]);
% 	update = dH(1:end-1);

	% restricted mean curvature flow.
% 	update = zeros(numel(H)-1,1,'single');
% 	sel = H(1:end-1) >= 0;
% 	update(sel) = max(dK1([sel; false]),0);
% 	sel = H(1:end-1) < 0;
% 	update(sel) = min(dK2([sel; false]),0);

	% S weighted dK1-dK2 flow.
% 	update = (dK2.*(S+1) + dK1.*(S-1)) / 2;
% 	update(end) = [];

	% relative movement weighted dK1-dK2 flow.
% 	update = zeros(size(K1),'single');
% 	sel = (abs(dK1) + abs(dK2)) > 0;
% 	update(~sel) = 0;
% 	update(sel)  = (dK1(sel).*abs(dK1(sel)) + dK2(sel).*abs(dK2(sel)))./(abs(dK1(sel)) + abs(dK2(sel)));
% 	update(end) = [];
	
	% restrict update to a maximum of some value so that things don't blow
	% up.
	update = max(min(update, 100), -100);

	update = update.*(denom.^0.5);

	change = max(abs(update));
	if( change < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	if( verbose >= 1 )
		fprintf('\n%4d: %6.4f, %6.4f, %f', iter, change, sqrt(mean(update.^2)), dt );
	end

	% Compute dt adaptively with each iteration.
	% Find the largest value of dt such that oscillations do not set in.
	oldphi = phi(i000);
% 	while(dt >= 0.0625)
		phi_new  = phi(i000) + dt*update;
% 		% strictly speaking I want to check for boundedness for +ve central
% 		% update only with neighboring max_phi resulting from decreasing
% 		% negative update and vice-versa.
% %  		if( all( min(backdrop(allNeighbors),[],2) <= phi_new & phi_new <= max(backdrop(allNeighbors),[],2)) )
% % 		if( all( (update < 0 & min(phi_more(allNeighbors), [], 2) <= phi_new) ...
% % 			   | (update > 0 & max(phi_less(allNeighbors), [], 2) >= phi_new) ))
% 		backdrop(i000) = update;
% 		update_neigh = backdrop(allNeighbors);
% 		backdrop(i000) = phi_new;
% % 		backdrop(i000) = phi(i000) + dt*update.*(update<0) - 1000*(update > 0);
% 		[maxPhi, imaxPhi] = max(backdrop(allNeighbors), [], 2);
% 		% We check for these conditions only in bandwidth - 2 of the level
% 		% set because the level set at bandwidth do not have phis computed
% 		% for all their neighbors. Strictly speaking only 1.732 is needed
% 		% in case of a true eucledian distance transform.
% 		if( any(update > 0 & any(update_neigh < 0, 2) & maxPhi < phi(i000) & phi_new > phi(allNeighbors(sub2ind(size(allNeighbors), (1:size(allNeighbors,1))', imaxPhi))) & abs(phi0) <= bandwidth - 2) )
% 			dt = dt/2;
% 			continue;
% 		end
% 
% 		[minPhi, iminPhi] = min(backdrop(allNeighbors), [], 2);
% % 		backdrop(i000) = phi(i000) + dt*update.*(update>0) + 1000*(update < 0);
% 		if( any(update < 0 & any(update_neigh > 0, 2) & minPhi > phi(i000) & phi_new < phi(allNeighbors(sub2ind(size(allNeighbors), (1:size(allNeighbors,1))', iminPhi))) & abs(phi0) <= bandwidth - 2) )
% 			dt = dt/2;
% 			continue;
%  		end
%  		break;
% 	end

	phi(i000) = max(min(phi_new, phi0 + upperThreshold), phi0 - lowerThreshold);
% 	dt = dt*2;	% Choose a larger dt for the next iteration.

	if( mod(iter, steps) == 0 )
		dumpStatus;
	end
end
if( verbose >= 1 )
	fprintf('\n');
end

fprintf('Last maximum error %d\n', max(abs(update)) );
Kfull = zeros(size(phi));
Hfull = zeros(size(phi));
Kfull(i000) = K(1:end-1);
Hfull(i000) = H(1:end-1);

	function dumpStatus
		persistent p;
		if( verbose >= 2 )
			if( drawIsoSurface )
				set(0,'CurrentFigure',fig1);
				if( ~isempty(p) )
					try
						delete(p);
					catch err
					end
				end
% 				backdrop(i000) = phi(i000) == phi0 + upperThreshold | phi(i000) == phi0 - lowerThreshold;
% 				[faces,verts,colors] = isosurface(phi, 1e-100, backdrop); 
% 				p = patch('Vertices', verts, 'Faces', faces, ... 
% 						'FaceVertexCData', colors, ... 
% 						'FaceColor','interp', ... 
% 						'edgecolor', 'interp');
				[faces,verts] = isosurface(phi, 1e-100); 
				p = patch('Vertices', verts, 'Faces', faces, ... 
						'edgecolor', 'none', 'FaceColor', [0.5 0.8 1.0]);
				title(sprintf('Iteration %d', iter));
			end

			set(0,'CurrentFigure',fig2);
			colormap(hot);
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
% 			backdrop(i000) = max(spacing(1:end-1,:),[],2);
			backdrop(i000) = iNZ(1:end-1);
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
