function [phi, Kfull, Hfull] = antiAlias3D_dH2(img, verbose, movementThreshold, bandwidth, maxIterations, spacing)
%
% Anti-aliases by using a dH flow.
%
% function [phi, Kfull, Hfull] = antiAlias3D_dH2(img)
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
end

if( nargin < 3 )
	% the absolute maximum movement allowed.
	% For true anti-aliasing, this should be 0.5, anything higher will
	% cause smoothing of data. Also, larger values of movementThreshold
	% should be accompanied by larger bandwidths.
	movementThreshold   = [0.5 0.5];
end

% the distance splat matrix.
persistent splat;
if( isempty(splat) )
	[xs, ys, zs] = ndgrid( -bandwidth-1:bandwidth+1, -bandwidth-1:bandwidth+1, -bandwidth-1:bandwidth+1 );
	splat = sqrt(xs.^2 + ys.^2 + zs.^2);
	clear xs ys zs
end

dt  = 0.02;
eps	= 0.01;			% if max(curvature) is less than this,
					% then assume that the level set
					% evolution has converged.

% Are we already passed in a distance map?
if( strcmp(class(img),'double') == 1 || strcmp(class(img),'single') == 1 )
	phi = single(img);
	inside = img < -0.5;
	if( verbose >= 1 )
		fprintf('Using passed-in distance map\n');
	end
else
	if( ~islogical(img) )
		img = logical(img);
	end
	inside = img;
	
	% TODO: Speed this up using our own distancing function.
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
lastPhiAtRedistance = phi0;

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

ordinates = [x y z];
clear x y z;

% the threshold for every vertex
upperThreshold = zeros(numel(i000), 1, 'single');
lowerThreshold = zeros(numel(i000), 1, 'single');

% build neighbor list for every vertex.
ineighbors = zeros(numel(i000)+1, 6, 'uint32');
neighCount = zeros(numel(i000)+1, 1, 'single');
iAllNeighbors = zeros(numel(i000), 26, 'uint32');
allNeighCount = zeros(numel(i000), 1, 'single');
% create a hash table to do a quick reverse lookup in i000.
ji000 = HashMap(numel(i000), 0.5);
for i = 1:numel(i000)
	ji000.put(i000(i), uint32(i));
end
% A temporary variable to store neighbors (make sure this is a row vector).
neighs = zeros(1,26, 'uint32');
for vi = 1:numel(i000)
	% 6, 18, or 26 neighbor mask -- choose.
 	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi)];
  	otherPts = [ ...
 		imm0(vi), ipp0(vi), imp0(vi), ipm0(vi), ...
 		im0m(vi), ip0p(vi), im0p(vi), ip0m(vi), ...
 		i0mm(vi), i0pp(vi), i0mp(vi), i0pm(vi), ...
		immm(vi), immp(vi), impm(vi), impp(vi), ipmm(vi), ipmp(vi), ippm(vi), ippp(vi) ];
	upperThreshold(vi) = movementThreshold(1);  % inside
	lowerThreshold(vi) = movementThreshold(2);  % outside  
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
	ineighbors(vi,:) = neighs(1:6);
	neighCount(vi) = neighsCount;
	% now continue on for all neighbors
	for ipt = 1:numel(otherPts)
		% do a reverse lookup for the index of the point.
		ni = ji000.get(otherPts(ipt));
		if( ~isempty(ni) )
			neighsCount = neighsCount + 1;
			neighs(neighsCount) = ni;
		end
	end
	neighs(neighsCount+1:end) = numel(i000)+1;
	iAllNeighbors(vi,:) = neighs;
	allNeighCount(vi) = neighsCount;
end
ineighbors(end,:) = numel(i000)+1;
clear neighs neighsCount vi pt ipt ni immm immp impm impp ipmm ipmp ippm ippp;

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
 	xd = ceil(size(phi,2)/2);
	yd = 1:size(phi,1);
 	zd = 1:size(phi,3);

% 	xd = 420;
% 	yd = 198:320;
% 	zd = 29:53;

% 	xd = 299;
% 	yd = 300:500;
% 	zd = 70:90;

	initZeroLevelSet = contour(squeeze(phi(yd,xd,zd)), 1e-100);
	initZeroLevelSet = initZeroLevelSet(:,2:end);
end
steps = 11;

% mean curvature - pre-allocated.
% the one extra element is to store a zero.
H = zeros(numel(i000)+1,1,'single');
K = H;
iNZ = false(size(H));

spacing = ones([numel(H) 3]);
distances = ones(size(ineighbors));
update = zeros(size(i000));

% the freezing mask. Prevent thin regions from evolving.
mask = true(size(i000));
% To keep track of the current state of a voxel. The check is done only if
% the voxel's crosses the actual value of the current level it is in.
lastInsideLevel = 100;

for iter = 1:maxIterations
%	if( any(abs(phi(i000) - lastPhiAtRedistance) >= 0.8) )
%		redistance;
%	end
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

	discriminant = H.^2 - 2*K;
	% discriminant can sometimes become slightly less than zero due to
	% numerical errors.
	discriminant = max(discriminant,0);
% 	sel = discriminant >= 0;
% 	assert( all(sel) );
% 	K1 = zeros(size(H), 'single');
% 	K2 = K1;
% 	K1(sel) = (H(sel) + sqrt(discriminant(sel)))*0.5;
% 	K2(sel) = (H(sel) - sqrt(discriminant(sel)))*0.5;
	K1 = (H + sqrt(discriminant))*0.5;
	K2 = (H - sqrt(discriminant))*0.5;

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
% 	update = zeros(numel(H)-1,1,'single');
% 	sel = K(1:end-1) < 0;
% 	update(sel) = dH([sel; false]);
% 	sel = K(1:end-1) >= 0 & H(1:end-1) > 0;
%  	update(sel) = max(dK1([sel; false]),0);
% % 	update(sel) = dK1([sel; false]);
% 	sel = K(1:end-1) >= 0 & H(1:end-1) < 0;
%  	update(sel) = min(dK2([sel; false]),0);
% % 	update(sel) = dK2([sel; false]);

	% Relaxed mean curvature flow.
% 	update = dH(1:end-1);

	% restricted mean curvature flow. K2 or K1 is selected depending upon
	% whether K1 is greater or K2 is greater.
% 	update = zeros(numel(H)-1,1,'single');
% 	sel = H(1:end-1) >= 0;
% 	update(sel) = dK1([sel; false]);
% 	sel = H(1:end-1) < 0;
% 	update(sel) = dK2([sel; false]);

	% alternate dK1 - dK2 flow
	if( mod(iter,2) == 0 )
		update = dK1(1:end-1);
	else
		update = dK2(1:end-1);
	end
	
 	% restricted mean curvature flow -- convex points do not flow outwards
 	% and concave points do not flow inwards.
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
	
	% convert update to a velocity.
	update = update.*(denom.^0.5);

	% find out areas to be frozen -- only update masks for newly added
	% voxels. One should be a tad more liberal with updating masks but it
	% would slow things down considerably.
	insideLevel = phi0 <= floor(phi0 + 0.5);
	for voxel = find(lastInsideLevel ~= insideLevel)'
		% TODO: add code to look for all neighbors and check if they need
		% to be frozen.
		for ni = 1:allNeighCount(voxel)
			npt = iAllNeighbors(voxel,ni);
			x = ordinates(npt,1);
			y = ordinates(npt,2);
			z = ordinates(npt,3);
			mask(npt) = ~toBeFrozen( phi(x-1:x+1,y-1:y+1,z-1:z+1) <= phi(i000(npt)) + 0.5 );
		end
	end
	lastInsideLevel = insideLevel;

	% for frozen areas, evolution is only stalled in the inward direction.
	% inward evolution happens by positive updates to the distance map.
	update(~mask) = min(0, update(~mask));
	old_phi = phi(i000);
  	phi_new = max(min(old_phi + dt*update, phi0 + upperThreshold), phi0 - lowerThreshold);
% 	phi_new = old_phi + dt*update;
	phi(i000) = phi_new;
	change  = max(abs((phi_new - old_phi)))/dt;

	if( change < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	if( verbose >= 1 )
		fprintf('\n%4d: %6.4f, %6.4f, %f', iter, change, sqrt(mean(update.^2)), dt );
	end
	
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

	function redistance
		%
		% redistances the level set.
		%
		redistance = false;
		if( verbose >= 1 )
			fprintf(' ... redistancing ... ' );
		end
		in  = i000( phi(i000) < -sqrt(3)/2 );
		out = i000( phi(i000) > sqrt(3)/2 );
		% update distances
		%phi(in)  = -inf;
		%phi(out) = inf;
		phi(in)  = -1000; % some large value instead of infinity
		phi(out) = 1000; % some large value instead of infinity
		% update the inside/outside mask.
		inside(in)  = true;
		inside(out) = false;
		level0  = abs(phi(i000)) <= sqrt(3)/2;	% the zero level set.
		% redistance the level set.
		for ilevel0pt = find(level0)'
			% a padding of bandwidth around the image would make things so much
			% simpler here :(
			% the entire x,y,z range.
			xs = [ordinates(ilevel0pt,1) - bandwidth-1, ordinates(ilevel0pt,1) + bandwidth+1];
			ys = [ordinates(ilevel0pt,2) - bandwidth-1, ordinates(ilevel0pt,2) + bandwidth+1];
			zs = [ordinates(ilevel0pt,3) - bandwidth-1, ordinates(ilevel0pt,3) + bandwidth+1];
			% the clipped x,y,z range.
			xc = [ max(xs(1), 1), min(xs(2),sz(1)) ];
			yc = [ max(ys(1), 1), min(ys(2),sz(2)) ];
			zc = [ max(zs(1), 1), min(zs(2),sz(3)) ];
			% "deposit" distances, negative distances inside and positive outside.
			modifier = inside( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ) * (-2) + 1;
			phi( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ) = modifier .* ...
				min( modifier .* phi( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ), ...
				splat( (xc(1)-xs(1)+1):(xc(2)-xs(1)+1), (yc(1)-ys(1)+1):(yc(2)-ys(1)+1), (zc(1)-zs(1)+1):(zc(2)-zs(1)+1) ) ...
				+ modifier .* phi( ordinates(ilevel0pt,1), ordinates(ilevel0pt,2), ordinates(ilevel0pt,3) ) );
		end
		if( verbose >= 1 )
			fprintf(' done.' );
		end
		lastPhiAtRedistance = phi(i000);
		% Ideally, we should also recompute i000 based on bandwidth,
		% however, doing so will be very expensive.
	end

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
% 				backdrop(i000) = phi(i000) == phi0 + upperThreshold | phi(i000) == phi0 - lowerThreshold | ~mask;
% 				backdrop(i000) = change;
% 				[faces,verts,colors] = isosurface(phi, 1e-100, backdrop); 
% 				p = patch('Vertices', verts, 'Faces', faces, ... 
% 						'FaceVertexCData', colors, ... 
% 						'FaceColor','interp', ... 
% 						'edgecolor', 'interp');
				[faces,verts] = isosurface(phi, 1e-100); 
				p = patch('Vertices', verts, 'Faces', faces, ... 
						'edgecolor', 'none', 'FaceColor', [0.5 0.8 1.0], ...
						'FaceLighting', 'gouraud');
				title(sprintf('Iteration %d', iter));
			end

			set(0,'CurrentFigure',fig2);
			colormap(hot);
			clf;

			subplot(2,2,1);
			backdrop(i000) = H(1:end-1);
			imagesc(squeeze(backdrop(yd,xd,zd)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			c = contour(squeeze(phi(yd,xd,zd)), 1e-100,'k');
			c = c(:,2:end);

			subplot(2,2,2);
% 			backdrop(i000) = max(spacing(1:end-1,:),[],2);
			backdrop(i000) = mask;
			imagesc(squeeze(backdrop(yd,xd,zd)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,3);
			backdrop(i000) = K(1:end-1);
			imagesc(squeeze(backdrop(yd,xd,zd)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,4);
			backdrop(i000) = (phi_new - old_phi)/dt;
			imagesc(squeeze(backdrop(yd,xd,zd)), [-change change]);
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			drawnow;
		end
	end
end


function [freeze] = toBeFrozen(inFG)

persistent freezeTable;

if( isempty(freezeTable) )
	if( exist('freezeTable.mat', 'file') )
		load freezeTable;
	else
		freezeTable = false(2^26,1);
		% A padded foreground is used to prevent neighbor out-of-bounds
		% checking.
		fg = true(5,5,5);
		nonPaddedFg = true(3,3,3);
		% foreground is true, background is false.
		powers_of_2 = pow2(1-26:0);
		one = uint32(1);
		% A temporary stack to contain the region while we grow it.
		voxels = zeros(27);
		% offsets for the 6 neighbors in the padded array.
		% being conservative and checking for only 6 neighbors instead of
		% all 26 ensures that development is frozen for even thin diagonal
		% regions.
		offsets = [ -1 +1 -5 +5 -25 +25 ];
		fprintf('Computing freezeLookupTable %6.2f%%', 0.0);
 		for i = 0:2^26-1
			if( mod(i, 10000) == 0 )
				fprintf('\b\b\b\b\b\b\b%6.2f%%', i/2^26*100);
			end
			mask = bitand( uint32(i * powers_of_2 - 0.5), one );
			nonPaddedFg([1:13 15:27]) = mask;
			fg(2:end-1,2:end-1,2:end-1) = nonPaddedFg;
			% start with 1st bg pixel.
			pos = find(~fg,1);
			if( isempty(pos) )
				% hmm, a hole will be formed if this voxel switches,
				% shouldn't happen.
				freeze = true;
			else
				% the central pixel will never be included because it is
				% in the foreground.
				voxels(1) = pos;
				voxelsTop = 1;
				while( voxelsTop > 0 )
					% pop a voxel, find its bg neighbors, push them on the
					% stack, reset this voxel to visited (by making it fg)
					pos = voxels(voxelsTop);
					voxelsTop = voxelsTop - 1;
					fg(pos) = true;
					neighbors = pos + offsets;
					% no need to check for out-of-bound neighbors because
					% of the padding.
					pos = neighbors(~fg(neighbors));
					newVoxelsTop = voxelsTop + numel(pos);
					voxels(voxelsTop+1:newVoxelsTop) = pos;
					voxelsTop = newVoxelsTop;
				end
				% if any bg pixels are left, then it means that there are
				% two disjoint bg regions. Development should be frozen for
				% such a region.
				freeze = any(~fg(:));
			end
			freezeTable(i+1) = freeze;
		end
		fprintf('\b\b\b\b\b\b\b%6.2f%%\n', 100);
		save freezeTable.mat freezeTable;
	end
end

% lookup in the freezeTable
mask = inFG([1:13 15:27]);
pow2_26 = pow2(25:-1:0);
freeze = freezeTable(sum(mask .* pow2_26)+1);

end
