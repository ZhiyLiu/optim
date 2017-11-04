function [phi, K, H] = antiAlias3D_dH(img, verbose, movementThreshold, maxIterations, spacing)
%
% Anti-aliases by using a dH flow.
%
% function phi = antiAlias3D_dH(img, verbose, movementThreshold, maxIterations, spacing)
%	Anti aliases a 3D binary image (1 inside, 0 outside and the boundary
% assumed to lie at 0.5); return a level-set with 0 representing the
% boundary.
%
% Input:
%	img					The input binary image (a 3D array)
%	verbose				0 for no progress, 1 for text, 2 for text+visuals
%						(default).
%	movementThreshold	0.5 for anti-aliasing, more for smoothing, a 2D
%						array, the first for inside, the second for outside
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
if( nargin < 3 )
	% the absolute maximum movement allowed.
	% For true anti-aliasing, this should be 0.5, anything higher will
	% cause smoothing of data. Also, larger values of movementThreshold
	% should be accompanied by larger bandwidths.
	movementThreshold   = 0.5;
end
spread    = 2;  % the size of the stencils (2*spread+1)
% this is the area in which curvatures will be computed
bandwidth = movementThreshold + sqrt(3)/2 + sqrt(2);
% the distance splat matrix.
persistent splat;
if( isempty(splat) )
	[xs, ys, zs] = ndgrid( -spread-1:spread+1, -spread-1:spread+1, -spread-1:spread+1 );
	splat = sqrt(xs.^2 + ys.^2 + zs.^2);
	clear xs ys zs
end

if( nargin < 4 )
	maxIterations	= 5000;
end
if( nargin < 5 )
	spacing = [1 1 1];
else
	spacing = spacing ./ min(spacing);
end

dt  = 0.02;
eps	= 0.001;		% if max(velocity) is less than this,
					% then assume that the level set
					% evolution has converged.

% Are we already passed in a distance map?
if( strcmp(class(img),'double') == 1 || strcmp(class(img),'single') == 1 )
	phi    = single(img);
	inside = img < -0.5;
	if( verbose >= 1 )
		fprintf('Using passed-in distance map\n');
	end
else
	if( ~islogical(img) )
		img    = logical(img);
	end
	inside = img;
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
% An image to store mean curvatures.
H = zeros(size(phi), 'single');

% Find the mask for the band of the computation.
img = -bandwidth(1) <= phi & phi <= bandwidth(2);
% mask off all boundary regions
img(1:spread, :, :) = false;
img(:, 1:spread, :) = false;
img(:, :, 1:spread) = false;
img(end-spread+1:end, :, :) = false;
img(:, end-spread+1:end, :) = false;
img(:, :, end-spread+1:end) = false;

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

ordinates = [x y z];
clear x y z;

% build neighbor list for every vertex.
ineighbors = zeros(numel(i000), 26, 'uint32');
neighCount = zeros(numel(i000), 1, 'single');
% create a hash table to do a quick reverse lookup in i000.
ji000 = java.util.HashMap(numel(i000), 0.5);
for i = 1:numel(i000)
	ji000.put(i000(i), uint32(i));
end
% A temporary variable to store neighbors (make sure this is a row vector).
neighs = zeros(1, 26, 'uint32');
for vi = 1:numel(i000)
  	pt = [im00(vi), ip00(vi), i0m0(vi), i0p0(vi), i00m(vi), i00p(vi) ...
 		imm0(vi), ipp0(vi), imp0(vi), ipm0(vi), ...
 		im0m(vi), ip0p(vi), im0p(vi), ip0m(vi), ...
 		i0mm(vi), i0pp(vi), i0mp(vi), i0pm(vi), ...
		immm(vi), immp(vi), impm(vi), impp(vi), ipmm(vi), ipmp(vi), ippm(vi), ippp(vi) ];
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
clear neighs neighsCount vi pt ipt ni ji000;

% the threshold for every point
upperThreshold = repmat(single(movementThreshold(1)), [numel(i000), 1]);	% inside
lowerThreshold = repmat(single(movementThreshold(2)), [numel(i000), 2]);	% outside

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
%  	xd = ceil(size(phi,2)/2);
% 	yd = 1:size(phi,1);
%  	zd = 1:size(phi,3);
%  	xd = 1:size(phi,2);
% 	yd = 1:size(phi,1);
%  	zd = ceil(size(phi,3)/2);

	xd = 420;
	yd = 198:320;
	zd = 29:53;

% 	xd = 299;
% 	yd = 300:500;
% 	zd = 70:90;

	initZeroLevelSet = contour(squeeze(phi(yd,xd,zd)), 1e-100);
	initZeroLevelSet = initZeroLevelSet(:,2:end);
end
steps = 21;

% spacing
spacing = [ 1 1 1 ];
spacing = repmat(spacing,[numel(i000) 1]);
redistance = false;

% the freezing mask. Prevent thin regions from evolving.
mask = true(size(i000));
lastLevel0 = false(size(i000));
lastPhiAtRedistance = phi0;
lastInStatus = phi(i000) <= 0;

for iter = 1:maxIterations
	% find pixels adjoining the 0 level set.
	level0  = abs(phi(i000)) <= sqrt(3)/2;	% zone to be updated
	level1  = abs(phi(i000)) <= sqrt(3)/2 + sqrt(2);	% zone where curvatures will be computed.
	update  = zeros(size(level0));
	%iNZ = false([numel(level1) 1]);
	% set the redistance flag, if the level set has moved more than some
	% amount less than a voxel. Do not use gradient to measure this as
	% the gradient can get close to zero at shocks.
	if( any(abs(phi(i000) - lastPhiAtRedistance) >= 0.8) )
		redistance = true;
	end
	% pixels that are not on the zero level set and inside or outside.
	if( redistance )
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
		% redistance the level set.
		for ipt = find(level0)'
			% a padding of spread around the image would make things so much
			% simpler here :(
			% the entire x,y,z range.
			xs = [ordinates(ipt,1) - spread-1, ordinates(ipt,1) + spread+1];
			ys = [ordinates(ipt,2) - spread-1, ordinates(ipt,2) + spread+1];
			zs = [ordinates(ipt,3) - spread-1, ordinates(ipt,3) + spread+1];
			% the clipped x,y,z range.
			xc = [ max(xs(1), 1), min(xs(2),sz(1)) ];
			yc = [ max(ys(1), 1), min(ys(2),sz(2)) ];
			zc = [ max(zs(1), 1), min(zs(2),sz(3)) ];
			% "deposit" distances, negative distances inside and positive outside.
			modifier = inside( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ) * (-2) + 1;
			phi( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ) = modifier .* ...
				min( modifier .* phi( xc(1):xc(2), yc(1):yc(2), zc(1):zc(2) ), ...
				splat( (xc(1)-xs(1)+1):(xc(2)-xs(1)+1), (yc(1)-ys(1)+1):(yc(2)-ys(1)+1), (zc(1)-zs(1)+1):(zc(2)-zs(1)+1) ) ...
				+ modifier .* phi( ordinates(ipt,1), ordinates(ipt,2), ordinates(ipt,3) ) );
		end
		if( verbose >= 1 )
			fprintf(' done.' );
		end
		lastPhiAtRedistance = phi(i000);
	end

	% now compute phi derivatives
	phix	= (phi(ip00(level1)) - phi(im00(level1))) ./ ( 2*spacing(level1,1) );
	phiy	= (phi(i0p0(level1)) - phi(i0m0(level1))) ./ ( 2*spacing(level1,2) );
	phiz	= (phi(i00p(level1)) - phi(i00m(level1))) ./ ( 2*spacing(level1,3) );
	phixx	= (phi(im00(level1)) - 2*phi(i000(level1)) + phi(ip00(level1))) ./ (spacing(level1,1).^2);
	phiyy	= (phi(i0m0(level1)) - 2*phi(i000(level1)) + phi(i0p0(level1))) ./ (spacing(level1,2).^2);
	phizz	= (phi(i00m(level1)) - 2*phi(i000(level1)) + phi(i00p(level1))) ./ (spacing(level1,3).^3);
	phixy	= (phi(ipp0(level1)) + phi(imm0(level1)) - phi(ipm0(level1)) - phi(imp0(level1))) ./ (4*spacing(level1,1).*spacing(level1,2));
	phiyz	= (phi(i0pp(level1)) + phi(i0mm(level1)) - phi(i0pm(level1)) - phi(i0mp(level1))) ./ (4*spacing(level1,2).*spacing(level1,3));
	phixz	= (phi(ip0p(level1)) + phi(im0m(level1)) - phi(ip0m(level1)) - phi(im0p(level1))) ./ (4*spacing(level1,1).*spacing(level1,3));

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom   = phix_2 + phiy_2 + phiz_2;

	% compute curvatures:
	% if we are redistancing all the time, then the gradient should always
	% approximately be equal to 1 (except near shocks)
	% non zero elements (greater then some small positive number)
	iNZ = denom > 0.1;
	level1 = find(level1);
	% mean curvature
	H(i000(level1(~iNZ))) = 0;
	H(i000(level1(iNZ)))  = ...
		 ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
		+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
		+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
		- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
		- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
		- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

	% compute curvature derivatives.
	Hx	= (H(ip00(level0)) - H(im00(level0))) ./ ( 2*spacing(level0,1) );
	Hy	= (H(i0p0(level0)) - H(i0m0(level0))) ./ ( 2*spacing(level0,2) );
	Hz	= (H(i00p(level0)) - H(i00m(level0))) ./ ( 2*spacing(level0,3) );
	Hxx	= (H(im00(level0)) - 2*H(i000(level0)) + H(ip00(level0))) ./ (spacing(level0,1).^2);
	Hyy	= (H(i0m0(level0)) - 2*H(i000(level0)) + H(i0p0(level0))) ./ (spacing(level0,2).^2);
	Hzz	= (H(i00m(level0)) - 2*H(i000(level0)) + H(i00p(level0))) ./ (spacing(level0,3).^3);
	Hxy	= (H(ipp0(level0)) + H(imm0(level0)) - H(ipm0(level0)) - H(imp0(level0))) ./ (4*spacing(level0,1).*spacing(level0,2));
	Hyz	= (H(i0pp(level0)) + H(i0mm(level0)) - H(i0pm(level0)) - H(i0mp(level0))) ./ (4*spacing(level0,2).*spacing(level0,3));
	Hxz	= (H(ip0p(level0)) + H(im0m(level0)) - H(ip0m(level0)) - H(im0p(level0))) ./ (4*spacing(level0,1).*spacing(level0,3));

	% non-zero in level 1 array
	l0in1 = level0(level1);
	% non-zero in level 0 array
	iNZ1  = denom > 0.1 & l0in1;
	iNZ   = denom(l0in1) > 0.1;

	% apply the laplacian-beltrami on H.
	dH = zeros( sum(level0), 1, 'single' );
	dH(iNZ) = ((Hyy(iNZ) + Hzz(iNZ)).*phix_2(iNZ1) ... 
			 + (Hxx(iNZ) + Hzz(iNZ)).*phiy_2(iNZ1) ...
			 + (Hxx(iNZ) + Hyy(iNZ)).*phiz_2(iNZ1) ...
			 - 2*phix(iNZ1).*phiy(iNZ1).*Hxy(iNZ) ...
			 - 2*phix(iNZ1).*phiz(iNZ1).*Hxz(iNZ) ...
			 - 2*phiy(iNZ1).*phiz(iNZ1).*Hyz(iNZ)) ./ denom(iNZ1) ...
			- (Hx(iNZ).*phix(iNZ1) + Hy(iNZ).*phiy(iNZ1) + Hz(iNZ).*phiz(iNZ1)) ...
			 .*H(iNZ) ./ denom(iNZ1).^0.5;

	update = -dH;
	% restrict update to a maximum of some value so that things don't blow
	% up
	update = max(min(update, 100), -100);

	% multiply by gradient of distance map.
	update = update.*(denom(l0in1).^0.5);

	% find out areas to be frozen -- only update masks for newly added
	% voxels. One should be a tad more liberal with updating masks but it
	% would slow things down considerably.
	inStatus = phi(i000) <= 0.0;
	for voxel = find(level0 & ~lastLevel0 | lastInStatus ~= inStatus )'
		% TODO: add code to look for all neighbors and check if they need
		% to be frozen.
		for ni = 1:neighCount(voxel)
			x = ordinates(ineighbors(voxel,ni),1);
			y = ordinates(ineighbors(voxel,ni),2);
			z = ordinates(ineighbors(voxel,ni),3);
			mask(ineighbors(voxel,ni)) = ~toBeFrozen( phi(x-1:x+1,y-1:y+1,z-1:z+1) <= 0 );
		end
	end
	lastInStatus = inStatus;
	
	% for frozen areas, evolution is only stalled in the inward direction.
	% inward evolution happens by positive updates to the distance map.
	% level0 voxels in mask -- mask(level0)
	update(~mask(level0)) = min(0, update(~mask(level0)));
	lastLevel0 = level0;
	old_phi = phi(i000(level0));
  	phi_new = max(min(old_phi + dt*update, phi0(level0) + upperThreshold(level0)), phi0(level0) - lowerThreshold(level0));
% 	phi_new = old_phi + dt*update;
	
	phi(i000(level0)) = phi_new;
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


H(i000(~level1)) = 0;
K = H;
K(i000(level1(iNZ)))  = ...
	 (phix_2(iNZ).*(phiyy(iNZ).*phizz(iNZ) - phiyz(iNZ).^2) ...
	+ phiy_2(iNZ).*(phixx(iNZ).*phizz(iNZ) - phixz(iNZ).^2) ...
	+ phiz_2(iNZ).*(phixx(iNZ).*phiyy(iNZ) - phixy(iNZ).^2) ...
	+ 2*(phix(iNZ).*phiy(iNZ).*(phixz(iNZ).*phiyz(iNZ) - phixy(iNZ).*phizz(iNZ)) ...
	   + phiy(iNZ).*phiz(iNZ).*(phixy(iNZ).*phixz(iNZ) - phiyz(iNZ).*phixx(iNZ)) ...
	   + phix(iNZ).*phiz(iNZ).*(phixy(iNZ).*phiyz(iNZ) - phixz(iNZ).*phiyy(iNZ)))) ...
	./ denom(iNZ).^2;


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
 				backdrop(i000(~level0)) = 0;
 				backdrop(i000(level0)) = (phi_new-old_phi)/dt;
				[faces,verts,colors] = isosurface(phi, 1e-100, backdrop); 
				p = patch('Vertices', verts, 'Faces', faces, ... 
						'FaceVertexCData', colors, ... 
						'FaceColor','interp', ... 
						'edgecolor', 'interp');
% 				[faces,verts] = isosurface(phi, 1e-100); 
% 				p = patch('Vertices', verts, 'Faces', faces, ... 
% 						'edgecolor', 'none', 'FaceColor', [0.5 0.8 1.0], ...
% 						'FaceLighting', 'gouraud');
				title(sprintf('Iteration %d', iter));
			end

			set(0,'CurrentFigure',fig2);
			colormap(hot);
			clf;
			title(sprintf('Iteration %d', iter));

			subplot(2,2,1);
			imagesc(double(squeeze(H(yd,xd,zd))));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			c = contour(squeeze(phi(yd,xd,zd)), 1e-100,'k');
			c = c(:,2:end);

			subplot(2,2,2);
% 			backdrop(i000(~level0))= 0;
% 			backdrop(i000(level0)) = phi_new >= phi0(level0) + upperThreshold(level0) | phi_new <= phi0(level0) - lowerThreshold(level0);
			backdrop(i000) = mask;
			imagesc(squeeze(backdrop(yd,xd,zd)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,3);
			backdrop(i000) = phi(i000);
			imagesc(squeeze(backdrop(yd,xd,zd)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			plot(c(1,:), c(2,:), 'k');

			subplot(2,2,4);
			backdrop(i000(~level0))= 0;
			backdrop(i000(level0)) = (phi_new-old_phi)/dt;
			range = min(change, 1);
			imagesc(squeeze(backdrop(yd,xd,zd)), [-range range]);
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
