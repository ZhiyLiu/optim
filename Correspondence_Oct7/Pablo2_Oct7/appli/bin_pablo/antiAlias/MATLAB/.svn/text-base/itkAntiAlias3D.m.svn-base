function phi = itkAntiAlias3D(img, verbose)
%
% function phi = itkAntiAlias3D(img)
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
maxIterations		= 1000;
bandwidth			= 6;

% .. and other not so parameter-like constants.
dt					= 0.0625;		% should be less than (1/2^(#dims+1))
movementThreshold   = 0.5;
eps = 0.05;							% if max(curvature) is less than this,
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
	fig1 = figure;
	view([-37.5 30]);
	fig2 = figure;
end

% the intial level set eroded and dilated by half a voxel.
phi0m05 = phi(i000) - 0.5;
phi0p05 = phi(i000) + 0.5;

if( verbose >= 1 )
	fprintf('Iterations     ');
end
for iter = 0:maxIterations
	if( verbose >= 1 )
		fprintf('\b\b\b\b%4d', iter);
	end
	if( mod(iter, (maxIterations/100)) == 0 && verbose >= 2)
		set(0,'CurrentFigure',fig1);
		[az, el] = view;
		clf;
		isosurface(phi,1e-100);
		axis vis3d;
		axis on;
		box on;
		axis equal;
		view(az,el);

		set(0,'CurrentFigure',fig2);
		clf;
		z = ceil(size(phi,3)/2);
		imagesc(phi(:,:,z), [-4 4]);
		axis image;
		colorbar;
		hold on;
		contour(phi(:,:,z),1e-100,'k');

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
	% non zero elements.
	iNZ = denom > 0;
	% mean curvature
	H(iNZ) = ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
			+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
			+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
			- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
			- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
			- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

	if( max(abs(H)) < eps )
		fprintf('\nLevel set is not changing.');
		break;
	end
	phi(i000)  = min(max(phi(i000) + dt*H, phi0m05), phi0p05);
end
if( verbose >= 1 )
	fprintf('\n');
end

fprintf('Last maximum error %d\n', max(abs(H)) );
