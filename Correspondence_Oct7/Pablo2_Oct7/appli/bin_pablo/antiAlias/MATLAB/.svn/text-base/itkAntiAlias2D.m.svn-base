function phi = itkAntiAlias2D(img)
% inside ~= 0, outside == 0
img = logical(img);
dout = bwdist( img ); din = -bwdist( ~img );
% correct for the fact that we assume pixel values to be at the center of
% the cell and the contour that is to be extracted is the one for 0.5
% (between 0 and 1)
dout = dout - 0.5;
din  = din + 0.5;
% put them together to one distance map
dall = zeros( size( img ) );
dall( img ) = din( img );
dall( ~img ) = dout( ~img );

nrIter = 5000;
dt = 0.1;

% Initialize,

phi = dall; z = zeros( size( img ) );
phix = z; phiy = z; phixx = z; phiyy = z; phixy = z;
figure(1);
figure(2);

dx = 1;
dy = 1;

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
		axis equal; box on;
		title(sprintf('Iteration %4d', iter));

		% draw the contour overlaid on the image.
		set(0,'CurrentFigure',2);
		clf
		subplot(1,2,1);
 		imshow(img, [0,1]);
		hold on;
		contour( phi, [0 0], 'r', 'linewidth', 1 );
		axis equal; box on;
		title('Original Image');

		subplot(1,2,2);
 		imshow(phi, [-5,5]);
		hold on;
		contour( phi,[0 0], 'r', 'linewidth', 1 );
		axis equal; box on;
		title('img');
		drawnow;
	end
	
% 	phi  = min(max(phi + dt*update, dall-0.5),dall+0.5);
 	phi  = phi + dt*update;
end

end