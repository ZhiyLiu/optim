function phi = antiAlias3D(img, verbose, movementThreshold, bandwidth, maxIterations, spacing)
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

import java.util.HashMap;

if( nargin < 2 )
	verbose = 2;
end

% some parameters
if( nargin < 4 )
	bandwidth		= 4;
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
% .. and other not so parameter-like constants.
dt					= 0.0625;		% should be less than (1/2^(#dims+1))
eps = 0.1;							% if max(curvature) is less than this,
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

% We will use a separate array to disable computations for those regions
% that have reached the stopping criteria for evolution to avoid
% recomputing the indices array -- evolveMask.

% find the indices for the mask, and then find the indices for all the
% derivatives. The x and y usage here is consistent with ndgrid instead of
% meshgrid, but that doesn't matter. Make sure it is a row vector;
i000 = uint32(find(img))';

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

clear x y z sz;

% build neighbor list for every vertex.
neighbors = cell(size(i000));

% create a hash table to do a quick reverse lookup in i000.
ji000 = HashMap(numel(i000), 0.5);
for i = 1:numel(i000)
	ji000.put(i000(i), uint32(i));
end
% A temporary variable to store neighbors (make sure this is a row vecotr).
% change to size 18 if using the 18 neighbor mask.
neighs = zeros(1,6, 'uint32');
for vi = 1:numel(i000)
	% 6 or 18 neighbor mask?
	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi)];
%  	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi) ...
% 		imm0(vi), ipp0(vi), imp0(vi), ipm0(vi), ...
% 		im0m(vi), ip0p(vi), im0p(vi), ip0m(vi), ...
% 		i0mm(vi), i0pp(vi), i0mp(vi), i0pm(vi) ];
	neighsCount = 0;
	for ipt = 1:numel(pt)
		% do a reverse lookup for the index of the point.
		ni = ji000.get(pt(ipt));
		if( ~isempty(ni) )
			neighsCount = neighsCount + 1;
			neighs(neighsCount) = ni;
		end
	end
	neighbors{vi} = neighs(1:neighsCount);
end
clear ji000 neighs neighsCount vi pt ipt ni;

if( verbose >= 1 )
	fprintf('Computed band and its indices.\n');
end

if( verbose >= 2 )
	fig1 = figure(1);
	set(0,'CurrentFigure',fig1);
	clf;
	origSurface = patch(isosurface(phi, 1e-100));
	set(origSurface, 'FaceColor', [0 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
	view([331 20]);
	camlight; lighting gouraud;
	xlabel('x');
	ylabel('y');
	zlabel('z');
	axis vis3d;
	axis on;
	box on;
	daspect(1./spacing);

	fig2 = figure(2);
	backdrop = zeros(size(phi),'single');
	
	% Added these variables here to facilitate debugging.
	drawIsoSurface = true;
	steps = 1;
% 	xs = ceil(size(phi,2)/2);
% 	ys = 1:size(phi,1);
% 	zs = 1:size(phi,3);
	xs = 420;
	ys = 198:320;
	zs = 29:53;
end

% the intial level set.
phi0 = phi(i000);


% A mask to prevent convex and concave regions from over-evolving
evolveMask = true(size(i000));

% dt changes from region to region to ensure that each region can evolve
% upto it's full potential and is not stuck because of a large choice of
% the initial dt.
adaptiveDt = repmat(dt, [1 numel(phi0)]);


% a temporary stack
stack = zeros(size(i000),'uint32');

regionId = zeros(size(i000),'uint32');
regions = struct('vertices', 1, 'neighbors', 1, 'freeRegionPtr', 1);
regions.vertices  = cell(size(i000));
regions.neighbors = cell(size(i000));
% build the free list pointer.
for ri = 1:numel(i000)-1
	regions.vertices{ri} = ri+1;
end
regions.vertices{numel(i000)} = 0;

% a temporary array to keep track of region as it is growing.
region = zeros(size(i000), 'uint32');

spacing = [1 1 1];	% reset spacing

iter = 0;
while iter < maxIterations
	iter = iter + 1;
	phix	= (phi(ip00) - phi(im00)) / ( 2*spacing(1) );
	phiy	= (phi(i0p0) - phi(i0m0)) / ( 2*spacing(2) );
	phiz	= (phi(i00p) - phi(i00m)) / ( 2*spacing(3) );
	phixx	= (phi(im00) - 2*phi(i000) + phi(ip00)) / (spacing(1).^2);
	phiyy	= (phi(i0m0) - 2*phi(i000) + phi(i0p0)) / (spacing(2).^2);
	phizz	= (phi(i00m) - 2*phi(i000) + phi(i00p)) / (spacing(3).^3);
	phixy	= (phi(ipp0) + phi(imm0) - phi(ipm0) - phi(imp0)) / (4*spacing(1)*spacing(2));
	phiyz	= (phi(i0pp) + phi(i0mm) - phi(i0pm) - phi(i0mp)) / (4*spacing(2)*spacing(3));
	phixz	= (phi(ip0p) + phi(im0m) - phi(ip0m) - phi(im0p)) / (4*spacing(1)*spacing(3));

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom	= phix_2 + phiy_2 + phiz_2;

	H = zeros(size(denom),'single');
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

	discriminant = H.^2 - 2*K;
	sel = discriminant > 0;
	K1 = zeros(size(H), 'single');
	K2 = K1;
	K1(sel) = (H(sel) + sqrt(discriminant(sel)))*0.5;
	K2(sel) = (H(sel) - sqrt(discriminant(sel)))*0.5;

	% classify all points
	ptClass = zeros(size(H),'uint8');	% unassigned region.
	ptClass( H < 0 & K > 0 ) = 1;		% concave
	ptClass( K <= 0 )		 = 2;		% hyperbolic + parabolic + flat regions.
	ptClass( H > 0 & K > 0 ) = 3;		% convex
	
	if( all(K >= 0 & H >= 0) || all(K >= 0 & H <= 0) )
        % This is an extremely rare condition to reach. It prevents an
        % ellipsoid from collapsing.
		fprintf('\nLevel set has been convexified.');
		break;
	end

	% if any vertex changes it's ptClass, then nullify the region and any
	% adjacent regions with this or the old ptClass.
	if( exist('lastPtClass', 'var') )
		numRegions = 0;
		for vi = find(lastPtClass ~= ptClass)
			ri = regionId(vi);
			% TODO: add a heuristic to decide if region needs to  be
			% nullified or not.
			% ok so the problem is that some other points could have been
			% changed, so we can't really just can't unassign everything
			% and free regions here.
			numRegions = numRegions + 1;
			region(numRegions) = ri;
			% check neighbors to see if they have the same class.
			for ni = neighbors{vi}
				ri = regionId(ni);
				% check if region has already been nullified.
				if( ri ~= 0 && (ptClass(vi) == ptClass(ni) || ...
						lastPtClass(vi) == ptClass(ni) || ...
						lastPtClass(vi) == lastPtClass(ni)) )
					% nullify this region too.
					numRegions = numRegions + 1;
					region(numRegions) = ri;
				end
			end
		end
		for ri = unique(region(1:numRegions))
			regionId(regions.vertices{ri}) = 0;	% unassign all vertices.
			% free this region and include it in the free region list.
			regions.vertices{ri}  = regions.freeRegionPtr;
			regions.freeRegionPtr = ri;
			% remove this region from neighbor list of adjacent regions.
			for nri = regions.neighbors{ri}
				neighs = regions.neighbors{nri};
				neighs(neighs == ri) = [];
				regions.neighbors{nri} = neighs;
			end
			regions.neighbors{ri} = [];
		end
		
		% Reset time step for any point whose class has changed or that
		% which is now hyperbolic.
		adaptiveDt(lastPtClass ~= ptClass | ptClass == 2) = dt;
		% Also, reset the evolution mask.
		evolveMask(lastPtClass ~= ptClass) = true;
	end
	lastPtClass = ptClass;

	% find regions for all unassigned vertices.
	visited  = false( size(i000) );
	for vi = 1:numel(regionId)
		if( regionId(vi) == 0 )
			% get an unassigned region for this vertex.
			ri = regions.freeRegionPtr;
			regions.freeRegionPtr = regions.vertices{ri};
			regions.vertices{ri} = [];
			regions.neighbors{ri} = [];
			regionSize  = 0;

			visited(vi) = true;

			% A stack to do a depth-first search.
			npts        = 1;
			stack(npts) = vi;
			
			% neighbor region list.
			neighs = zeros(size(i000),'uint32');
			neighCount = 0;
			while( npts > 0 )
				pt   = stack(npts);
				npts = npts - 1;
				regionSize = regionSize + 1;
				region(regionSize) = pt;
				regionId(pt) = ri;
				% look at neighbors and add them if they are part of the
				% region.
				for i = neighbors{pt}
					% i is now index in i000 array for the pt.
					if( abs(phi(i000(i)) - phi(i000(vi))) <= 0.5 )
						if( ~visited(i) && ptClass(i) == ptClass(vi) && regionId(i) == 0 )
							% part of connected region, push the point in
							% the queue.
							visited(i)   = true;
							npts = npts + 1;
							stack(npts) = i;
						elseif( ptClass(i) ~= ptClass(vi) && regionId(i) ~= 0 )
							% add region corresponding to this vertex to
							% list of neighbors (and the other way too)
							% -- only if it has already been assigned a
							% regionId.
							nri = regionId(i);
							if( numel(regions.neighbors{nri}) == 0 ...
								|| regions.neighbors{nri}(end) ~= ri )
								% if we did not add this region to the
								% neighbor list of the neighbor, add it.
								regions.neighbors{nri}(end+1) = ri;
							end
							neighCount = neighCount + 1;
							neighs(neighCount) = nri;
						end
					end
				end
			end
			regions.vertices{ri}  = region(1:regionSize);
			if( neighCount > 0 )
				regions.neighbors{ri} = unique(neighs(1:neighCount));
			else
				regions.neighbors{ri} = [];
			end
		end
	end
	
	% now flow convex surfaces inwards by K2, concave surfaces outwards by
	% K1, hyperbolic surfaces connected only to convex surfaces outwards by
	% K2, connected only to concave surfaces inwards by K1, and all
	% remaining surfaces by mean curvature H (is H the best answer?).
	update = zeros(size(i000), 'single');
 	update(ptClass == 1) = K2(ptClass == 1);	% concave
	update(ptClass == 3) = K1(ptClass == 3);	% convex
	% hyperbolic and others.
	for vi = find(ptClass == 2)
		% Does this region have convex neighbors?
		convex  = false;
		% How about concave neighbors?
		concave = false;
		for nri = regions.neighbors{ regionId(vi) }
			if( ptClass( regions.vertices{nri}(1) ) == 1 )
				concave = true;
			elseif( ptClass( regions.vertices{nri}(1) ) == 3 )
				convex = true;
			else
				% ignore hyperbolic neighbors ...
			end
			if( convex && concave )
				break;
			end
		end
		if( convex && concave )
			update(vi) = H(vi);
		elseif( convex )
			% We just set an upper bound of -0.1 to speed things up a bit.
			update(vi) = min(K2(vi),-0.1);
		elseif( concave )
			% We just set a lower bound of 0.1 to speed things up a bit.
			update(vi) = max(K1(vi),0.1);
		end
		% there could be some hyperbolic regions with no neighbors, this is
		% not strictly correct in a geometric sense, but can happen because
		% of the way we construct regions.
	end
	update(~iNZ) = 0;
	update(iNZ) = update(iNZ) .* (denom(iNZ).^0.5);

	error = max(abs(update(evolveMask)));
	if( sum(evolveMask) == 0 || error < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	fprintf('\n%4d: %6.4f, %6.4f', iter, error, sqrt(mean(update(evolveMask).^2)) );
	
	if( verbose >= 2 && mod(iter, steps) == 0)
		if( drawIsoSurface )
			set(0,'CurrentFigure',fig1); 
			%[az, el] = view;
			if( exist('p', 'var') );
				delete(p);
			end
			p = patch(isosurface(phi,1e-100));
			set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
			%view(az,el);
		end

		set(0,'CurrentFigure',fig2);
		clf;

		subplot(2,2,1);
		backdrop(i000) = evolveMask + 1;
		imagesc(squeeze(backdrop(ys,xs,zs)), [0 2]);
% 		axis image;
		colorbar;
		hold on;
		contour(squeeze(phi(ys,xs,zs)),1e-100,'k');
		
		backdrop(i000) = H;
		subplot(2,2,2);
		imagesc(squeeze(backdrop(ys,xs,zs)));
% 		axis image;
		colorbar;
		hold on;
		contour(squeeze(phi(ys,xs,zs)),1e-100,'k');

		backdrop(i000) = K;
		subplot(2,2,3);
		imagesc(squeeze(backdrop(ys,xs,zs)));
% 		axis image;
		colorbar;
		hold on;
		contour(squeeze(phi(ys,xs,zs)),1e-100,'k');

		backdrop(i000) = ptClass;
		subplot(2,2,4);
		imagesc(squeeze(backdrop(ys,xs,zs)), [0 3]);
% 		axis image;
		colorbar;
		hold on;
		contour(squeeze(phi(ys,xs,zs)),1e-100,'k');

		drawnow;
% 		pause(0.1);
	end
	
	newphi = phi(i000) + adaptiveDt.*update;
	% Any pixel that has been displaced by more than a voxel (threshold)
	% should stop evolving. Any connecting region of this pixel,
	% connectivity given by K <= 0 will also stop evolving.
 	elements = find( (abs(newphi - phi0) > movementThreshold) & evolveMask);
	visited = false(size(i000));
	for ei = elements;
		if( ~visited(ei) )
			% find the region to which this pixel belongs.
			if( ptClass(ei) == 1 || ptClass(ei) == 3 )
				% do all this only for convex or concave regions.
				pts = regions.vertices{regionId(ei)};
			else
				pts  = pt;
			end
			visited(pts) = true;
			newDt = (movementThreshold - abs(phi0(pts) - phi(i000(pts)))) ./ abs(update(pts));
			if( any(newDt < 0.001) )
				dtRatio = 0;
			else
				% We multilpy by 0.9 to make the flow a little conservative
				% so that we don't overshoot.
				dtRatio = 0.9 * min(newDt./adaptiveDt(pts));
			end
			if( dtRatio == 0 )
				evolveMask(pts) = false;
			else
				adaptiveDt(pts) = dtRatio * adaptiveDt(pts);
			end
		end
	end
	phi(i000(evolveMask)) = phi(i000(evolveMask)) + adaptiveDt(evolveMask).*update(evolveMask);
end
if( verbose >= 1 )
	fprintf('\n');
end

fprintf('Last maximum error %d\n', max(abs(update)) );
