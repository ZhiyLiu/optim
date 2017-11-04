function phi = antiAlias3D_reflectedFlow(img, verbose, movementThreshold, bandwidth, maxIterations, spacing)
%
% Uses alternate min(K2,0), max(K1,0) flow. The flow alternates only when
% any piece of the zero-level set hits a bound.
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

% find the indices for the mask, and then find the indices for all the
% derivatives. The x and y usage here is consistent with ndgrid instead of
% meshgrid, but that doesn't matter.
i000 = uint32(find(img));

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

if( verbose >= 1 )
	fprintf('Computed band and its indices.\n');
end

if( verbose >= 2 )
	fig1 = figure(1);
	set(0,'CurrentFigure',fig1);
	clf;
	origSurface = patch(isosurface(phi, 1e-100));
	set(origSurface, 'FaceColor', [0 1 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
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
	steps = 20;
% 	xs = ceil(size(phi,2)/2);
% 	ys = 1:size(phi,1);
% 	zs = 1:size(phi,3);
	xs = 420;
	ys = 198:320;
	zs = 29:53;

	initZeroLevelSet = contour(squeeze(phi(ys,xs,zs)), 1e-100);
	initZeroLevelSet = initZeroLevelSet(:,2:end);
end

% the intial level set.
phi0 = phi(i000);

flowDir = +1;
lastReversal = 1;
lastPass = +inf;
iter = 0;
spacing = [1 1 1];
while iter < maxIterations || any(iter - lastReversal == [floor(lastPass/2) ceil(lastPass/2)])
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
% 	K = H;
	% non zero elements.
	iNZ = denom > 0;
	% mean curvature
	H(iNZ) = ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
			+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
			+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
			- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
			- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
			- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

% 	% gaussian curvature
% 	K(iNZ) = (phix_2(iNZ).*(phiyy(iNZ).*phizz(iNZ) - phiyz(iNZ).^2) ...
% 			+ phiy_2(iNZ).*(phixx(iNZ).*phizz(iNZ) - phixz(iNZ).^2) ...
% 			+ phiz_2(iNZ).*(phixx(iNZ).*phiyy(iNZ) - phixy(iNZ).^2) ...
% 		+ 2*(phix(iNZ).*phiy(iNZ).*(phixz(iNZ).*phiyz(iNZ) - phixy(iNZ).*phizz(iNZ)) ...
% 		   + phiy(iNZ).*phiz(iNZ).*(phixy(iNZ).*phixz(iNZ) - phiyz(iNZ).*phixx(iNZ)) ...
% 		   + phix(iNZ).*phiz(iNZ).*(phixy(iNZ).*phiyz(iNZ) - phixz(iNZ).*phiyy(iNZ)))) ...
% 		./ denom(iNZ).^2;

% 	discriminant = H.^2 - 2*K;
% 	sel = discriminant > 0;
% 	K1 = zeros(size(H));
% 	K2 = K1;
% 	K1(sel) = (H(sel) + sqrt(discriminant(sel)))*0.5;
% 	K2(sel) = (H(sel) - sqrt(discriminant(sel)))*0.5;

% 	% alternate K1-K2 flow.
% 	if( flowDir > 0 )
% 		update = max(K1, 0.001);
% 	else
% 		update = min(K2, -0.001);
% 	end
% 	update = update .* (denom.^0.5);

	% alternate mean curvature flow
	update = flowDir * H.*(denom.^0.5);

	error = max(abs(update));
	if( error < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	if( verbose >= 1 )
		fprintf('\n%4d: %6.4f, %6.4f', iter, error, sqrt(mean(update.*update)) );
	end

%  	newphi = max(min(phi(i000) + dt*update, phi0 + 5), phi0 - 5);
 	newphi = max(min(phi(i000) + dt*update, phi0 + movementThreshold+1e-6), phi0 - movementThreshold-1e-6);
	clamped = any(abs(newphi) <= 1.0 & abs(newphi - phi0) > movementThreshold);
	if(clamped)
% 		% try with a smaller dt.
% 		sel = (movementThreshold+1e-6 - abs(phi0 - phi(i000))) ./ abs(update);
% 		newDt = min(sel(sel >= 0));
% 		assert(newDt < dt);
% 		if( newDt > 0.001 )
% 			if( verbose >= 1 )
% 				fprintf(' ... using a smaller dt of %d.', newDt);
% 			end
% 		 	phi(i000) = phi(i000) + newDt*update;
% 		else
			flowDir = -flowDir;
			lastPass = iter - lastReversal;
			if( verbose >= 1 )
				fprintf(' ... reversing flow (iterations %d)', lastPass);
			end
			dumpStatus;
			% Change i to see where it has clamped.
			% backdrop(i000) = (abs(newphi) <= 1.0 & abs(newphi - phi0) > movementThreshold) + 1; clamped = backdrop; backdrop(i000) = phi0;
			% for i=38; clf; imagesc(squeeze(clamped(ys,i,zs)),[0 2]); title(num2str(i)); hold on; contour(squeeze(backdrop(ys,i,zs)),1e-100,'w'); contour(squeeze(phi(ys,i,zs)),1e-100,'k', 'LineWidth', 2); pause(0.1); end
			lastReversal = iter;
% 		end
	else
		phi(i000) = newphi;
	end

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
				if( exist('p', 'var') );
					delete(p);
				end
				p = patch(isosurface(phi,1e-100));
				set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
			end

			set(0,'CurrentFigure',fig2);
			clf;
			backdrop(i000) = H;
			imagesc(squeeze(backdrop(ys,xs,zs)));
			colorbar;
			hold on;
			plot(initZeroLevelSet(1,:), initZeroLevelSet(2,:), 'w');
			contour(squeeze(phi(ys,xs,zs)),1e-100,'k', 'LineWidth', 2);

			drawnow;
		end
	end
end