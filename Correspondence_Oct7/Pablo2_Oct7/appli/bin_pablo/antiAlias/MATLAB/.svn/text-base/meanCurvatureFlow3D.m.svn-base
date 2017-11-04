function phi = meanCurvatureFlow3D(img, verbose)
%
% function phi = meanCurvatureFlow3D(img)
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
maxIterations		= 500;
% .. and other not so parameter-like constants.
dt					= 0.0625;		% should be less than (1/2^(#dims+1))

% compute the distance map.
dout = bwdist( img );
din = -bwdist( 1-img );
% correct for the fact that we assume pixel values to be at the center of
% the cell and the contour that is to be extracted is the one for 0.5
% (between 0 and 1) (matlab bwdist thinks differently)
% dout = dout - 0.5;
% din  = din + 0.5;
% put them together to one distance map
mask = logical(img);
phi = zeros( size( mask ), 'double' );
phix = phi;
phiy = phi;
phiz = phi;
phixx = phi;
phiyy = phi;
phizz = phi;
phixy = phi;
phiyz = phi;
phixz = phi;
phi( mask ) = din( mask );
phi( ~mask ) = dout( ~mask );
clear dout din
if( verbose >= 1 )
	fprintf('Computed initial level set.\n');
end

if( verbose >= 2 )
	fig1 = figure(1);
	view([-37.5 30]);
	fig2 = figure(2);
end

dx = 1; dy = 1; dz = 1;

% % Find the mask for the band of the computation.
% bandwidth = 5;
% mask = -bandwidth < phi & phi < bandwidth;
% % mask off all boundary regions
% mask([1 end], :, :) = false;
% mask(:, [1 end], :) = false;
% mask(:, :, [1 end]) = false;

for iter = 0:maxIterations
	phix(2:end-1,:,:)	= (phi(3:end,:,:) - phi(1:end-2,:,:)) / (2*dx);
	phiy(:,2:end-1,:)	= (phi(:,3:end,:) - phi(:,1:end-2,:)) / (2*dy);
	phiz(:,:,2:end-1)	= (phi(:,:,3:end) - phi(:,:,1:end-2)) / (2*dz);
	phixx(2:end-1,:,:)	= (phi(1:end-2,:,:) - 2*phi(2:end-1,:,:) + phi(3:end,:,:)) / (dx*dx);
	phiyy(:,2:end-1,:)	= (phi(:,1:end-2,:) - 2*phi(:,2:end-1,:) + phi(:,3:end,:)) / (dy*dy);
	phizz(:,:,2:end-1)	= (phi(:,:,1:end-2) - 2*phi(:,:,2:end-1) + phi(:,:,3:end)) / (dz*dz);
	phixy(2:end-1,2:end-1,:)	= (phi(3:end,3:end,:) + phi(1:end-2,1:end-2,:) - phi(3:end,1:end-2,:) - phi(1:end-2,3:end,:)) / (4*dx*dy);
	phiyz(:,2:end-1,2:end-1)	= (phi(:,3:end,3:end) + phi(:,1:end-2,1:end-2) - phi(:,3:end,1:end-2) - phi(:,1:end-2,3:end)) / (4*dy*dz);
	phixz(2:end-1,:,2:end-1)	= (phi(3:end,:,3:end) + phi(1:end-2,:,1:end-2) - phi(3:end,:,1:end-2) - phi(1:end-2,:,3:end)) / (4*dx*dz);

	phix_2	= phix.^2;
	phiy_2	= phiy.^2;
	phiz_2	= phiz.^2;
	denom	= phix_2 + phiy_2 + phiz_2;

	H = zeros(size(denom));
	% non zero elements.
% 	iNZ = denom > 0 & mask;
	iNZ = denom > 0;
	% mean curvature
	H(iNZ) = ((phiyy(iNZ) + phizz(iNZ)).*phix_2(iNZ) ... 
			+ (phixx(iNZ) + phizz(iNZ)).*phiy_2(iNZ) ...
			+ (phixx(iNZ) + phiyy(iNZ)).*phiz_2(iNZ) ...
			- 2*phix(iNZ).*phiy(iNZ).*phixy(iNZ) ...
			- 2*phix(iNZ).*phiz(iNZ).*phixz(iNZ) ...
			- 2*phiy(iNZ).*phiz(iNZ).*phiyz(iNZ)) ./ denom(iNZ).^1.5;

	if( mod(iter, 10) == 0 && verbose >= 2)
		set(0,'CurrentFigure',fig1); 
		[az, el] = view;
		clf;
		isosurface(phi,1e-100);
		xlabel('x');
		ylabel('y');
		zlabel('z');
		axis vis3d;
		axis on;
		box on;
		axis equal;
		view(az,el);

		set(0,'CurrentFigure',fig2);
		clf;
		x = ceil(size(phi,2)/2);

		imagesc(squeeze(H(:,x,:)));
		axis image;
		colorbar;
		hold on;
		contour(squeeze(img(:,x,:)),0.5,'w');
		contour(squeeze(phi(:,x,:)),1e-100,'k');

		drawnow;
	end
% 	fprintf('%d: maximum mean curvature %d\n', iter, max(abs(H(mask))) );
 	fprintf('%d: maximum mean curvature %d\n', iter, max(abs(H(:))) );
	% mean curvature flow
% 	phi(mask) = phi(mask) + dt*H(mask);
	phi = phi + dt*H.*(denom.^0.5);
end

if( verbose >= 1 )
	fprintf('\n');
end
