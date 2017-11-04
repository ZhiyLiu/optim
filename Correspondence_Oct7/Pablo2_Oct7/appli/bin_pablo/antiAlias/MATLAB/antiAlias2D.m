function phi = antiAlias2D(img, dx, dy)

nrIter = 5000;
dt = 0.1;
threshold = 0.5;
curvature_threshold = -0.0;

% inside ~= 0, outside == 0
img = logical(img);
dout = bwdist( img ); din = -bwdist( ~img );
% correct for the fact that we assume pixel values to be at the center of
% the cell and the contour that is to be extracted is the one for 0.5
% (between 0 and 1)
dout = dout - 0.5;
din  = din + 0.5;
% put them together to one distance map
phi0 = zeros( size( img ) );
phi0( img ) = din( img );
phi0( ~img ) = dout( ~img );

% Initialize,
mask = true( size(img) );
newphi = phi0;
phi = phi0; z = zeros( size( img ) );
phix = z; phiy = z; phixx = z; phiyy = z; phixy = z;
figure(1);
figure(2);

for iter=0:nrIter
	% compute the derivatives
	phix(2:end-1,:) = (phi(3:end,:)-phi(1:end-2,:))/(2*dx);
	phiy(:,2:end-1) = (phi(:,3:end)-phi(:,1:end-2))/(2*dy);
	phixx(2:end-1,:) = (phi(1:end-2,:)-2*phi(2:end-1,:)+phi(3:end,:))/(dx*dx);
	phiyy(:,2:end-1) = (phi(:,1:end-2)-2*phi(:,2:end-1)+phi(:,3:end))/(dy*dy);
	phixy(2:end-1,2:end-1) = (phi(3:end,3:end)+phi(1:end-2,1:end-2)-phi(3:end,1:end-2)- ...
		phi(1:end-2,3:end))/(4*dx*dy);
	
	denom = phix.^2 + phiy.^2;
	iNZ = find( denom>0 );
	% When the denom is zero, the numerator will also be zero.
	update = z;
	update(iNZ) = (phixx(iNZ).*phiy(iNZ).^2-2*phix(iNZ).*phiy(iNZ).*phixy(iNZ)+ ...
		phiyy(iNZ).*phix(iNZ).^2)./denom(iNZ);

	if ( mod(iter,100)==0 )
		fprintf('%d: maximum curvature %d\n', iter, max(abs(update(:))) );
		set(0,'CurrentFigure',1);
		clf;
		contour( phi,1e-100, 'k', 'linewidth', 1 );
		box on;
		daspect( [dx dy 1] );
		title(sprintf('Iteration %4d', iter));

		% draw the contour overlaid on the image.
		set(0,'CurrentFigure',2);
		clf
		subplot(2,1,1);
 		imshow(img, [0,1]);
		hold on;
		contour( phi, 1e-100, 'r', 'linewidth', 1 );
		box on;
		daspect( [dx dy 1] );
		title('Original Image');

		subplot(2,1,2);
 		imshow(phi, [-5,5]);
		hold on;
		contour( phi, 1e-100, 'r', 'linewidth', 1 );
		box on;
		daspect( [dx dy 1] );
		title('img');
		drawnow;
	end
	
	newphi(mask) = phi(mask) + dt*update(mask);
	% Any pixel that has been displaced by more than a voxel should stop
	% evolving. Any connecting region of this pixel, connectivity given by
	% K <= 0 will also stop evolving.
	elements = find((abs(newphi - phi0) > threshold) & mask )';
	for i = elements
		[y,x] = ind2sub(size(mask),i);
		% if this pixel has not already been masked out by another point in
		% this loop, do a flood fill.
		if( mask(y,x) )
			% term 'abs(phi0-newphi(y,x)) <= 2' is to ensure that we don't
			% spill out into too far a level.
			mask(bwselect(update < curvature_threshold & mask & abs(phi0-newphi(y,x)) <= 0.5, x, y, 4)) = false;
			% ensure that this point becomes false in any case.
			mask(y,x) = false;
		end
	end
	phi(mask) = newphi(mask);

	% pure mean curvature flow
% 	phi  = min(max(phi + dt*update, phi0-0.5),phi0+0.5);
%  	phi  = phi + dt*update;
end

end

% load shapes
% get spiral shape from an image
% load tstSmall
% im = imread('aliasedRectangle.png');
% im = imread('weirdShape.png');
% im = imread('NonConvexPolygon.png');
% im = imread('DogBone.png');

